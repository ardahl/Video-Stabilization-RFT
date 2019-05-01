%  Video Stabilization using RANSAC method
%   Additions to the simple method are tested here before added to base.
function StabilizeVideoSimple2(video, output_name, varargin)
    dispMethodDefault = 0;
    p = inputParser;
    validPosNum = @(x) isnumeric(X) && isscalar(x) && x>=0;
    addOptional(p, 'dispMethod', dispMethodDefault, validPosNum);
    parse(p, varargin{:});

    vid = VideoReader(video);
    
    if p.Results.dispMethod == 0
        hVPlayer = vision.VideoPlayer;
    end
    numFrames = 0;
    fprintf('Counting Frames\n');
    while hasFrame(vid)
        readFrame(vid);
        numFrames = numFrames + 1;
        if numFrames ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%d', numFrames);
    end
    fprintf('\n%d Frames\n', numFrames);
    frames = cell(numFrames,1);
    framesGray = cell(numFrames,1);
    transforms = repmat(eye(3), 3*numFrames, 1);
    Tparams = repmat([1 0 0 0], numFrames, 1);
    vid.CurrentTime = 0;
    frame = readFrame(vid);
    frameg = rgb2gray(frame);
    frames{1} = frame;
    framesGray{1} = frameg;
    framewarp = frame;
    d1 = detectSURFFeatures(frameg);
    [f1, valid1] = extractFeatures(frameg, d1);
    ii = 2;
    while hasFrame(vid)
        img = frame;
        imgwarp = framewarp;
        frame = readFrame(vid);
        frameg = rgb2gray(frame);
        frames{ii} = frame;
        framesGray{ii} = frameg;

        % Collect Features from both images
        d2 = detectSURFFeatures(frameg);
        [f2, valid2] = extractFeatures(frameg, d2);
        % Match features between the 2 images
        indexPairs = matchFeatures(f1, f2);
        pts1  = valid1(indexPairs(:,1));
        pts2 = valid2(indexPairs(:,2));
        % Estimate transform using MSAC
        [H, inlier2, inlier1] = estimateGeometricTransform(pts2, pts1, 'affine');
        imgBp = imwarp(frame, H, 'OutputView', imref2d(size(frame)));
        pointsBmp = transformPointsForward(H, inlier2.Location);
        
        % Refit H as scale-rotation-translation for better numerical
        % stability
        [scale, theta, translation] = decomposeT(H.T);
        % Reconstitute new s-R-t transform:
        HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; translation], [0 0 1]'];
        tformsRT = affine2d(HsRt);
        
%         transforms(3*ii-2:3*ii,:) = HsRt;
        Tparams(ii,:) = [scale, theta, translation];
        
        imgBold = imwarp(frame, H, 'OutputView', imref2d(size(frame)));
        imgBsRt = imwarp(frame, tformsRT, 'OutputView', imref2d(size(frame)));
        
        framewarp = imwarp(frame,affine2d(HsRt),'OutputView',imref2d(size(frame)));
        
        if p.Results.dispMethod == 0
            step(hVPlayer, imfuse(imgwarp,framewarp,'ColorChannels','red-cyan'));
        elseif p.Results.dispMethod == 1
            figure;
            showMatchedFeatures(img,frame,pts1,pts2);
            waitforbuttonpress;
            close;
        elseif p.Results.dispMethod == 2
            figure;
            showMatchedFeatures(img, imgBp, inlier1, pointsBmp);
            waitforbuttonpress;
            close;
        elseif p.Results.dispMethod == 3
            figure;
            imshowpair(imgBold,imgBsRt,'ColorChannels','red-cyan');
            waitforbuttonpress;
            close;
        end
                
        f1 = f2;
        valid1 = valid2;
        ii = ii+1;
    end
    if p.Results.dispMethod == 0
        hide(hVPlayer);
        release(hVPlayer);
    end
    
    % Smooth the transformations
    fprintf('Smoothing\n');
    smoothTparams = Tparams;
    k = 3;
    for i=1:numFrames
        lb = max(1, i-k);
        ub = min(numFrames, i+k);
        neighbors = Tparams(lb:ub,:);
        ave = mean(neighbors);
        smoothTparams(i,:) = ave;
        scale = ave(1);
        theta = ave(2);
        translation = ave(3:4);
        transforms(3*i-2:3*i,:) = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; translation], [0 0 1]'];
    end
    
    sx = [-1 -2 -1; 0 0 0; 1 2 1];
    sy = [-1 0 1; -2 0 2; -1 0 1];
    blurWeights = zeros(numFrames, 2*k+1);
    bt = zeros(numFrames,1);
    % Compute relative blurriness for each image
    fprintf('Blurriness Metric\n');
    for i=1:numFrames
        if i ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*i/numFrames, i, numFrames);
        f = framesGray{i};
        f = im2double(f);
        b = sum(sum(conv2(f,sx).^2 + conv2(f,sy).^2));
        bt(i) = 1./b;
    end
    vid.CurrentTime = 0;
    % Compute weighting between image and neighbors
    fprintf('\nBlur Weighting\n');
    alpha = 1;
    for i=1:numFrames
        if i ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*i/numFrames, i, numFrames);
        b = bt(i);
        im1 = im2double(framesGray{i});
        lb = max(1, i-k);
        ub = min(numFrames, i+k);
        for j=lb:ub
            bp = bt(j);
            ratio = b / bp;
            if ratio > 1
                im2 = framesGray{j};
                im2 = im2double(imwarp(im2,affine2d(transforms(j*3-2:j*3,:)),'OutputView',imref2d(size(im2))));
                err = norm(im2-im1);
                blurWeights(i,j) = ratio * (alpha ./ (err + alpha));
            end
        end
    end
    
%     fprintf('Mosaicing\n');
%     fprintf('Computing Variance: ');
%     white = ones(size(framesGray{1}));
%     for ind=1:numFrames
%         if ind ~= 1
%             fprintf(repmat('\b', 1, nchar));
%         end
%         nchar = fprintf('%0.2f%% (%d/%d)', 100*ind/numFrames, ind, numFrames);
%         frame = frames{ind};
%         T = transforms(ind*3-2:ind*3,:);
%         emptyMask = imwarp(white,affine2d(T),'OutputView',imref2d(size(frame)));
%         emptyMask = imbinarize(emptyMask);
%         if ~any(any(emptyMask))
%             continue
%         end
%         lb = max(1, ind-k);
%         ub = min(numFrames, ind+k);
%         nNeigh = ub-lb+1;
%         % first calculate Ibar
%         aveIntensity = zeros(size(framesGray));
%         for i=lb:ub
%             img = frame{j};
%             Tn = transforms(i*3-2:i*3,:);
%             neighwarp = imwarp(img,affine2d(Tn),'OutputView',imref2d(size(img)));
%             intensity = 0.3*neighwarp(:,:,1) + 0.59*neighwarp(:,:,2) + 0.11*neighwarp(:,:,3);
%             aveIntensity = aveIntensity + intensity;
%         end
%         aveIntensity = aveIntensity ./ nNeigh;
%         % calculate variance
%         var = zeros(size(framesGray));
%         for i=lb:ub
%             Tn = transforms(i*3-2:i*3,:);
%             neighwarp = imwarp(frame,affine2d(Tn),'OutputView',imref2d(size(frame)));
%             intensity = 0.3*neighwarp(:,:,1) + 0.59*neighwarp(:,:,2) + 0.11*neighwarp(:,:,3);
%             diff = intensity - aveIntensity;
%             var = var + diff.^2;
%         end
%         var = var / (nNeigh-1);
%         Minds = find(emptyMask==true);
%         variance = var(Minds); %#ok<FNDSB>
%         for i=1:numel(variance)
%             
%         end
%     end
    
    fprintf('\nWriting Video\n');
    stable = VideoWriter(output_name);
    open(stable);
    maxDeblur = 0;
    maxDeblurNorm = 0;
    maxWeight = 0;
    for ind=1:numFrames
        if ind ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*ind/numFrames, ind, numFrames);
        frame = frames{ind};
        T = transforms(ind*3-2:ind*3,:);
        framewarp = imwarp(frame,affine2d(T),'OutputView',imref2d(size(frame)));
        framewarp = im2double(framewarp);
        weightSum = 0;
        pixelSum = zeros(size(frame));
        for i=ind-k:ind+k
            if i < 1 || i > numFrames
                continue
            end
            w = blurWeights(ind,i);
            weightSum = weightSum + w;
            Tp = transforms(i*3-2:i*3,:);
            framewarp2 = imwarp(frames{i},affine2d(Tp),'OutputView',imref2d(size(frames{i})));
            framewarp2 = w * im2double(framewarp2);
            pixelSum = pixelSum + framewarp2;
        end
        normtmp = norm(rgb2gray(pixelSum));
        frame_deblur = (framewarp + pixelSum) / (1 + weightSum);
        if weightSum > maxWeight
            maxDeblur = ind;
            maxDeblurNorm = normtmp;
            maxWeight = weightSum;
            before = framewarp;
            after = frame_deblur;
            mask = pixelSum;
        end
        writeVideo(stable, frame_deblur);
    end
    fprintf('\nDone\n');
    
    fprintf('Most deblurred frame: %f, %f, %f\n', maxDeblur, maxDeblurNorm, maxWeight);
    mon = cat(4, before, after);
    figure;
    montage(mon);
    figure;
    imshowpair(im2uint8(before),after,'falsecolor');
end

% Refit H as scale-rotation-translation for better numerical
        % stability
function [s, theta, translation] = decomposeT(T)
    R = T(1:2,1:2);
    % Compute theta from mean of two possible arctangents
    theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
    % Compute scale from mean of two stable mean calculations
    s = mean(R([1 4])/cos(theta));
    % Translation remains the same:
    translation = T(3, 1:2);
end