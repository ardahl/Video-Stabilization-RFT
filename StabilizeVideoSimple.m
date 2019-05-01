%  Video Stabilization using RANSAC method
% video is the name of the file to stabilize
% output_name is the name of the file to write the video to (will be .avi)
% Options
% DispMethod - 0 for live player showing the difference between current and
%       previous frame, 1/2 show matching features between different warps
%       of the frames, 3 for showing the difference between full affine and
%       similarity matrix
% Smooth - true if video should be smoothed
% SmoothK - How many neighbors on both side to use for average, total of
%       2*k+1 values
% ShowError - Plots the norm of the difference between averaging methods
%       and prints bounds of parameters
% ShowPath - Plots the x translation values with base and both smoothing
%       methods
function StabilizeVideoSimple(video, output_name, varargin)
    dispMethodDefault = 0;
    smoothDefault = true;
    smoothKDefault = 5;
    showErrorDefault = false;
    showPathDefault = false;
    p = inputParser;
    validPosNum = @(x) isnumeric(x) && isscalar(x) && x>=0;
    validLogical = @(x) isscalar(x) && islogical(x);
    addOptional(p, 'DispMethod', dispMethodDefault, validPosNum);
    addOptional(p, 'Smooth', smoothDefault, validLogical);
    addOptional(p, 'SmoothK', smoothKDefault, validPosNum);
    addOptional(p, 'ShowError', showErrorDefault, validLogical);
    addOptional(p, 'ShowPath', showPathDefault, validLogical);
    parse(p, varargin{:});

    vid = VideoReader(video);
    
    if p.Results.DispMethod == 0
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
    transforms = repmat(eye(3), 3*numFrames, 1);
    Tparams = repmat([1 0 0 0], numFrames, 1);
    vid.CurrentTime = 0;
    frame = readFrame(vid);
    frameg = rgb2gray(frame);
    framewarp = frame;
    d1 = detectSURFFeatures(frameg);
    [f1, valid1] = extractFeatures(frameg, d1);
    ii = 2;
    stable = VideoWriter(output_name);
    open(stable);
    while hasFrame(vid)
        img = frame;
        imgwarp = framewarp;
        frame = readFrame(vid);
        frameg = rgb2gray(frame);

        % Collect Features from both images
        d2 = detectSURFFeatures(frameg);
        [f2, valid2] = extractFeatures(frameg, d2);
        % Match features between the 2 images
        indexPairs = matchFeatures(f1, f2);
        pts1 = valid1(indexPairs(:,1));
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
        
        Tparams(ii,:) = [scale, theta, translation];
        
        imgBold = imwarp(frame, H, 'OutputView', imref2d(size(frame)));
        imgBsRt = imwarp(frame, tformsRT, 'OutputView', imref2d(size(frame)));
        
        framewarp = imwarp(frame,affine2d(HsRt),'OutputView',imref2d(size(frame)));
        
        if p.Results.DispMethod == 0
            step(hVPlayer, imfuse(imgwarp,framewarp,'ColorChannels','red-cyan'));
        elseif p.Results.DispMethod == 1
            figure;
            showMatchedFeatures(img,frame,pts1,pts2);
            waitforbuttonpress;
            close;
        elseif p.Results.DispMethod == 2
            figure;
            showMatchedFeatures(img, imgBp, inlier1, pointsBmp);
            waitforbuttonpress;
            close;
        elseif p.Results.DispMethod == 3
            figure;
            imshowpair(imgBold,imgBsRt,'ColorChannels','red-cyan');
            waitforbuttonpress;
            close;
        end
                
        f1 = f2;
        valid1 = valid2;
        ii = ii+1;
    end
    
    vid.CurrentTime = 0;
    ind = 1;
    if p.Results.Smooth || p.Results.ShowError
        k = p.Results.SmoothK;
    else
        k = 0;
    end
    if p.Results.ShowError
        norms = zeros(1,numFrames);
        sbounds = [Inf -Inf];
        tbounds = [Inf -Inf];
        xbounds = [Inf -Inf];
        ybounds = [Inf -Inf];
    end
    if p.Results.ShowPath
        path = zeros(numFrames,3);
    end
    fprintf('Writing Video\n');
    while hasFrame(vid)
        if ind ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*ind/numFrames, ind, numFrames);
        frame =  readFrame(vid);
        if p.Results.Smooth
            lb = max(1, ind-k);
            ub = min(numFrames, ind+k);
            neighbors = Tparams(lb:ub,:);
            ave = mean(neighbors,1);
            scale = ave(1);
            theta = ave(2);
            translation = ave(3:4);
            S = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; translation], [0 0 1]'];
            if p.Results.ShowPath
                path(ind,1) = Tparams(ind,3);
                path(ind,2) = translation(1);
                lb2 = max(1, ind-2*k);
                ub2 = min(numFrames, ind+2*k);
                neighbors2 = Tparams(lb2:ub2,:);
                ave2 = mean(neighbors2);
                path(ind,3) = ave2(3);
            end
            
            if p.Results.ShowError
                bnd = Tparams(ind,:);
                sbounds(1) = min(sbounds(1), bnd(1));
                sbounds(2) = max(sbounds(2), bnd(1));
                tbounds(1) = min(tbounds(1), bnd(2));
                tbounds(2) = max(tbounds(2), bnd(2));
                xbounds(1) = min(xbounds(1), bnd(3));
                xbounds(2) = max(xbounds(2), bnd(3));
                ybounds(1) = min(ybounds(1), bnd(4));
                ybounds(2) = max(ybounds(2), bnd(4));
                S2 = zeros(3,3);
                for n=lb:ub
                    vals = Tparams(n,:);
                    scale = vals(1);
                    theta = vals(2);
                    translation = vals(3:4);
                    Stmp = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; translation], [0 0 1]']; 
                    S2 = S2 + Stmp;
                end
                S2 = S2 ./ numel(lb:ub);
                
                norms(ind) = norm(S-S2, 'fro');
            end
        else
            ave = Tparams(ind,:);
            scale = ave(1);
            theta = ave(2);
            translation = ave(3:4);
            S = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; translation], [0 0 1]'];
        end
        framewarp = imwarp(frame,affine2d(S),'OutputView',imref2d(size(frame)));
        writeVideo(stable, framewarp);
        ind = ind + 1;
    end
    fprintf('\nDone\n');
    
    if p.Results.ShowPath
        figure
        hold on
        plot(path(:,1), 'g--', 'LineWidth', 1.5);
        plot(path(:,2), 'r-', 'LineWidth', 1.5);
        plot(path(:,3), 'b-', 'LineWidth', 1.5);
        title('X Translation');
        xlabel('frame');
        ylabel('x position');
        legend('Base',strcat('k=',num2str(k)),strcat('k=',num2str(2*k)));
    end
    
    if p.Results.ShowError
        plot(norms);
        title('Frobenious Norm of Average Components minus Average Matrices');
        xlabel('Frame');
        ylabel('Norm');
        fprintf('Max Error: %f\n', max(norms));
        fprintf('Unfiltered Bounds:\n');
        fprintf('scale: %f - %f\n', sbounds(1), sbounds(2));
        fprintf('theta: %f - %f\n', tbounds(1), tbounds(2));
        fprintf('dx: %f - %f\n', xbounds(1), xbounds(2));
        fprintf('dy: %f - %f\n', ybounds(1), ybounds(2));
    end
end

function [s, theta, translation] = decomposeT(T)
    R = T(1:2,1:2);
    % Compute theta from mean of two possible arctangents
    theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
    % Compute scale from mean of two stable mean calculations
    s = mean(R([1 4])/cos(theta));
    % Translation remains the same:
    translation = T(3, 1:2);
end