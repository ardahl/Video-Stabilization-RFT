% This code uses the Levenberg-Marquardt code by Alexander:
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox

%  Video Stabilization using Robust Feature Trajectories
% video is the name of the file to stabilize
% output_name is the name of the file to write the video to (will be .avi)
% Options
% SigmaWeighting - Sigma for the gaussian in the trajectory weights.
% SigmaTrajectory - (currently unused) Sigma for the gaussian in the trajectory cost
% MaxFrames - Only process the first n frames, -1 for the whole video
% MinTrajLength - All trajectories shorter than this will be removed
% LambdaT - Neighborhood consistancy constant in trajectory cost
% CostCutoff - Trajectories with cost above this will be retired
% ShowCosts - For debugging, shows graph of trajectory costs and exits
% ShowTrajectories - For debugging, graphs the path of all the trajectories
% FeatureType - Select either 'SURF', 'KAZE', or 'MSER' type features
function StabilizeVideoRFT(video, output_name, varargin)
    maxFramesDefault = -1;
    sigmaSDefault = 30;
    sigmaTDefault = 30;
    minTrajLengthDefault = 3;
    lambdaTDefault = 1;
    costCutoffDefault = 1000;
    showCostsDefault = false;
    showTrajDefault = false;
    p = inputParser;
    validPosNum = @(x) isnumeric(x) && isscalar(x) && x>=0;
    validNum = @(x) isnumeric(x) && isscalar(x);
    validLogical = @(x) isscalar(x) && islogical(x);
    validCutoff = @(x) validPosNum(x) || (ischar(x) && strcmp(x,'Outliers'));
    addOptional(p, 'SigmaWeighting', sigmaSDefault, validPosNum);
    addOptional(p, 'SigmaTrajectory', sigmaTDefault, validPosNum);
    addOptional(p, 'MaxFrames', maxFramesDefault, validNum);
    addOptional(p, 'MinTrajLength', minTrajLengthDefault, validPosNum);
    addOptional(p, 'LambdaT', lambdaTDefault, validPosNum);
    addOptional(p, 'CostCutoff', costCutoffDefault, validCutoff);
    addOptional(p, 'ShowCosts', showCostsDefault, validLogical);
    addOptional(p, 'ShowTrajectories', showTrajDefault, validLogical);
    addOptional(p, 'FeatureType', 'SURF', @(x)strcmp(x,'SURF')||strcmp(x,'KAZE')||strcmp(x,'MSER'));
    parse(p, varargin{:});

    vid = VideoReader(video);
    numFrames = 0;
    fprintf('Counting Frames\n');
    while hasFrame(vid)
        readFrame(vid);
        numFrames = numFrames + 1;
        if numFrames ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%d', numFrames);
        if numFrames == p.Results.MaxFrames
            break
        end
    end
    fprintf('\n%d Frames\n', numFrames);
    
    % Set of all trajectories
    % Each trajectory is a cell of 5 elements
        % [x, y] - a rowx2 vector containing the x,y path
        % [features] - row vector containing the feature descriptions at
            % each time
        % [ts, te] - a 1x2 vector containing the start and end frames
        % [neighbor indices] - array of indices of the neighboring traj.
        % live = int on whether the trajectory is considered live(1),
            % retired(2), or invalid(0)
    
    vid.CurrentTime = 0;
    
    % Process the first frame uniqely
    fprintf('\nComputing Feature Trajectories\n');
    nchar = fprintf('%0.2f%% (%d, %d)', 100*(1/numFrames), 1, numFrames);
    frame = readFrame(vid);
    frameg = rgb2gray(frame);
    frameg = im2double(frameg);
    sz = size(frameg);
    if strcmp(p.Results.FeatureType,'SURF')
        surfPoints = detectSURFFeatures(frameg);
        [features, validPoints] = extractFeatures(frameg, surfPoints);
    elseif strcmp(p.Results.FeatureType,'KAZE')
        kazePoints = detectKAZEFeatures(frameg);
        [features, validPoints] = extractFeatures(frameg, kazePoints, 'FeatureSize', 128);
    elseif strcmp(p.Results.FeatureType,'MSER')
        mserPoints = detectMSERFeatures(frameg);
        [features, validPoints] = extractFeatures(frameg, mserPoints);
    end
    fsize = size(features,2);
    % Extract x,y from features
    pts = double(validPoints.Location);
    % Allocate extra space so it doesn't resize as much
    % Assume average of length(features) features detected each frame
        % In order to reduce complexity at the expence of size, non-valid
        % trajectories of length 1 are kept.
    trajectories = cell(numFrames*length(features), 1);
    trajPts = zeros(numFrames*length(features), numFrames, 2);
    validPts = zeros(numFrames*length(features), numFrames);
    totalTrajUsed = length(features);
    initLength = numFrames*length(features);
    % Compute Delaunay triangulation from points
    tri = delaunayTriangulation(pts);
    edge = edges(tri);
    % Record neighbor information for each feature
    neigh = cell(size(pts,1),1);
    for i=1:size(edge,1)
        % Loop through each edge, need to mark down neighbor for (i,j) and
        % (j,i)
        e1 = edge(i,1);
        e2 = edge(i,2);
        neigh{e1} = [neigh{e1} e2];
        neigh{e2} = [neigh{e2} e1];
    end
    % Create new trajectories for each feature found
    for i=1:length(features)
        newT = cell(1,5);
        newT{1} = pts(i,:);
        trajPts(i,1,:) = pts(i,:);
        validPts(i,1) = 1;
        newT{2} = features(i,:);
        newT{3} = [1 1];
        newT{4} = neigh{i};
        newT{5} = 1;
        trajectories{i} = newT;
    end
    goodTraj = totalTrajUsed;
    
    pdTW = makedist('Normal',0,p.Results.SigmaTrajectory);
    if p.Results.ShowCosts
        costs = zeros(numFrames*length(features),2);
        numCosts = 0;
    end
    ind = 2;
    while hasFrame(vid)
        fprintf(repmat('\b', 1, nchar));
        nchar = fprintf('%0.2f%% (%d, %d)\r', 100*(ind/numFrames), ind, numFrames);
        
        frame = readFrame(vid);
        frameg = rgb2gray(frame);
        frameg = im2double(frameg);
        %% Addition
        % Collect Features from the new frame
        if strcmp(p.Results.FeatureType,'SURF')
            surfPoints = detectSURFFeatures(frameg);
            [features, validPoints] = extractFeatures(frameg, surfPoints);
        elseif strcmp(p.Results.FeatureType,'KAZE')
            kazePoints = detectKAZEFeatures(frameg);
            [features, validPoints] = extractFeatures(frameg, kazePoints, 'FeatureSize', 128);
        elseif strcmp(p.Results.FeatureType,'MSER')
            mserPoints = detectMSERFeatures(frameg);
            [features, validPoints] = extractFeatures(frameg, mserPoints);
        end
        %% Linking
        % Extract x,y from features
        pts = double(validPoints.Location);
        % Compute Delaunay triangulation from points
        tri = delaunayTriangulation(pts);
        edge = edges(tri);
        % Record neighbor information for each feature
        neigh = cell(size(pts,1),1);
        for i=1:size(edge,1)
            % Loop through each edge, need to mark down neighbor for (i,j) and
            % (j,i)
            e1 = edge(i,1);
            e2 = edge(i,2);
            neigh{e1} = [neigh{e1} e2];
            neigh{e2} = [neigh{e2} e1];
        end
%         trimesh(tri.ConnectivityList, tri.Points(:,1), tri.Points(:,2), ind+zeros(size(pts,1),1));
        
        %% Propigation
        % Match features between the 2 images
        % Recover the features of all live trajectories
        trajFeat = zeros(goodTraj,fsize,'single');
        trajInd = zeros(goodTraj);
        c1 = 1;
        for i=1:totalTrajUsed
            if isempty(trajectories{i})
                continue
            end
            if trajectories{i}{5} == 1
                trajFeat(c1,:) = trajectories{i}{2}(end,:);
                trajInd(c1) = i;
            end
            c1 = c1 + 1;
        end
        % IndexPairs(:,1) are the indices to the current trajectories 
        % (trajFeat) that are matching with the new features.
        % trajInd(indexPairs(:,1)) will give the indices to the trajectory
        % master list. They can then be updated with the indices to the new
        % features from indexPairs(:,2)
        indexPairs = matchFeatures(trajFeat, features, 'Unique', true);
        pts2 = pts(indexPairs(:,2),:);
        feat2 = features(indexPairs(:,2),:);
                
        % With the matching features, add them to the trajectories
        % Loop through each live trajectory and check if any of the
        % matching points from the previous from the previous frame are
        % part of the live set. 
        matchedTraj = trajInd(indexPairs(:,1));
        nodeInd = zeros(1,size(features,1));
        for i=1:size(pts2,1)
            T = trajectories{matchedTraj(i)};
            newPts = [T{1}; pts2(i,:)];
            trajPts(matchedTraj(i),ind,:) = pts2(i,:);
            validPts(matchedTraj(i),ind) = 1;
            newFeat = [T{2}; feat2(i,:)];
            trajectories{matchedTraj(i)}{1} = newPts;
            trajectories{matchedTraj(i)}{2} = newFeat;
            trajectories{matchedTraj(i)}{3}(2) = ind;
            nodeInd(indexPairs(i,2)) = matchedTraj(i);
        end
        
        % For the new matches, create a new trajectory
        feat2 = features;
        indMap = 1:size(features,1);
        feat2(indexPairs(:,2),:) = [];
        indMap(indexPairs(:,2)) = [];
        p2 = pts;
        p2(indexPairs(:,2),:) = [];
        for i=1:size(p2,1)
            newT = cell(1,5);
            newT{1} = p2(i,:);
            trajPts(totalTrajUsed+i,ind,:) = p2(i,:);
            validPts(totalTrajUsed+i,ind) = 1;
            newT{2} = feat2(i,:);
            newT{3} = [ind ind];
            newT{4} = [];
            newT{5} = 1;
            trajectories{totalTrajUsed+i} = newT;
            nodeInd(indMap(i)) = totalTrajUsed+i;
        end
        goodTraj = goodTraj + size(p2,1);
        totalTrajUsed = totalTrajUsed + size(p2,1);
        
        % Go through neighbor information and find neighboring trajectories
        % Loop through features2, find the corresponding feature in the
        % trajectory list. Build a list of indices where node index i gives 
        % the index to the trajectory.
        for i=1:size(features,1)
            % Find trajectory that this feature belongs to, mark node i as
            % part of that trajectory
            if nodeInd(i) < 1
                fprintf('I dont think this should happen\n');
                continue
            end
            nList = neigh{i};
            nt = zeros(numel(nList),1);
            for j=1:numel(nList)
                nt(j) = nodeInd(nList(j));
            end
            trajectories{nodeInd(i)}{4} = union(trajectories{nodeInd(i)}{4}, nt);
        end
        
        %% Pruning
        % Loop through current trajectories and find ones that didn't get
        % an addition, move to retired.
        for i=1:totalTrajUsed
            if isempty(trajectories{i})
                continue
            end
            if trajectories{i}{5} ~= 1
                continue
            end
            if trajectories{i}{3}(2) ~= ind && trajectories{i}{3}(2)-trajectories{i}{3}(1) < p.Results.MinTrajLength
                trajectories{i}{5} = 0; % length < n retired trajectories are invalid
                % Need at least length n to calculate acceleration for
                % optimization
                validPts(i,trajectories{i}{3}(2)) = 0;
            elseif trajectories{i}{3}(2) ~= ind
                trajectories{i}{5} = 2;
            end
        end
        
        % Calculate cost for each trajectory
        currCostInvalid = 0;
        if ischar(p.Results.CostCutoff)
            frameCosts = zeros(1,size(pts2,1));
        end
        for i=1:size(pts2,1)
            T = trajectories{matchedTraj(i)};
            siftDiff = norm(T{2}(size(T{2},1)-1,:)-T{2}(size(T{2},1),:)).^2;
            ui = T{1}(size(T{1},1),:) - T{1}(size(T{1},1)-1,:);
            nbor = T{4};
            nConsistancy = 0;
            for j=1:numel(nbor)
                nind = nbor(j);
                if isempty(trajectories{nind})
                    continue
                end
                Tj = trajectories{nind};
                if Tj{5} == 0 || Tj{3}(1) == Tj{3}(2)
                    continue
                end
                uj = Tj{1}(size(Tj{1},1),:) - Tj{1}(size(Tj{1},1)-1,:);
                flowNorm = norm(ui-uj).^2;
%                 lookBackDist = max(T{3}(1), Tj{3}(1));
%                 lookBackMax = min(ind-1, Tj{3}(2)-1);
%                 tau = numel(lookBackDist:ind-1);
%                 vNormSum = 0;
%                 for prevFrame=lookBackDist:lookBackMax
%                     vi = T{1}((prevFrame+1)-T{3}(1)+1,:) - T{1}(prevFrame-T{3}(1)+1,:);
%                     vj = Tj{1}((prevFrame+1)-Tj{3}(1)+1,:) - Tj{1}(prevFrame-Tj{3}(1)+1,:);
%                     vNormSum = vNormSum + norm(vi-vj).^2;
%                 end
%                 D = 1/tau * vNormSum;
%                 wij = pdf(pdTW, sqrt(D));
%                 nConsistancy = nConsistancy + wij*flowNorm;
                nConsistancy = nConsistancy + flowNorm;
            end
            cost = siftDiff + p.Results.LambdaT*nConsistancy;
            if p.Results.ShowCosts
                costs(numCosts+1,:) = [ind cost];
                numCosts = numCosts + 1;
            end
            % If the cost of the trajectory is above a threshold, prune it
            % and add it to the retired set
            if ~ischar(p.Results.CostCutoff) && cost > p.Results.CostCutoff
                trajectories{matchedTraj(i)}{5} = 2;
                trajectories{matchedTraj(i)}{1}(size(trajectories{matchedTraj(i)}{1},1),:) = [];
                trajectories{matchedTraj(i)}{2}(size(trajectories{matchedTraj(i)}{2},1),:) = [];
                trajectories{matchedTraj(i)}{3}(2) = ind - 1;
                if trajectories{matchedTraj(i)}{3}(2)-trajectories{matchedTraj(i)}{3}(1) < p.Results.MinTrajLength
                    trajectories{matchedTraj(i)}{5} = 0;
                    validPts(matchedTraj(i),trajectories{matchedTraj(i)}{3}(2)) = 0;
                end
                % Can ignore neighbor information as it won't be used for
                % anything anymore
                currCostInvalid = currCostInvalid + 1;
            elseif ischar(p.Results.CostCutoff)
                frameCosts(i) = cost;
            end
        end
        if ischar(p.Results.CostCutoff)
            TF = isoutlier(frameCosts,'quartiles','ThresholdFactor', 3);
            currCostInvalid = sum(TF);
            invalidInds = matchedTraj(TF);
            for i=1:numel(invalidInds)
                trajectories{invalidInds(i)}{5} = 2;
                trajectories{invalidInds(i)}{1}(size(trajectories{invalidInds(i)}{1},1),:) = [];
                trajectories{invalidInds(i)}{2}(size(trajectories{invalidInds(i)}{2},1),:) = [];
                trajectories{invalidInds(i)}{3}(2) = ind - 1;
                % Can ignore neighbor information as it won't be used for
                % anything anymore
                if trajectories{invalidInds(i)}{3}(2)-trajectories{invalidInds(i)}{3}(1) < p.Results.MinTrajLength
                    trajectories{invalidInds(i)}{5} = 0;
                    validPts(invalidInds(i),trajectories{invalidInds(i)}{3}(2)) = 0;
                end
            end
        end
        if size(pts2,1)-currCostInvalid < 10
            warning('Number of valid trajectories (%d) for frame %d is <10. This may result in lower quality tranformations. Try lowering the cost cutoff or the lambdaT parameters',size(pts2,1)-currCostInvalid,ind); 
        end
        
        % Clean up invalid ones
        numBad = 0;
        for i=1:totalTrajUsed
            if isempty(trajectories{i})
                continue
            end
            if trajectories{i}{5} == 0
                trajectories{i} = [];
                numBad = numBad + 1;
            end
        end
        goodTraj = goodTraj - numBad;
        
        ind = ind + 1;
        if ind == numFrames+1
            break
        end
    end
    
    if p.Results.ShowCosts
        if numCosts < initLength
            costs(numCosts+1:initLength,:) = [];
        end
        figure
        subplot(1,2,1);
        hold on
        % Get outlier points of the weights
        TF = isoutlier(costs(:,2),'quartiles', 'ThresholdFactor', 3);
        outliers = costs(TF,:);
        regular = costs(~TF,:);
        plot(regular(:,1),regular(:,2),'b.','LineWidth',1.5);
        plot(outliers(:,1),outliers(:,2),'r.','LineWidth',1.5);
        hold off
        
        subplot(1,2,2);
        hold on
        cind = 1;
        ol = [];
        reg = [];
        for i=1:numFrames
            secstart = cind;
            while cind <= size(costs,1) && costs(cind,1) == i
                cind = cind + 1;
            end
            section = costs(secstart:cind-1,:);
            TF = isoutlier(section(:,2),'quartiles', 'ThresholdFactor', 3);
            ol = [ol; section(TF,:)]; %#ok<AGROW>
            reg = [reg; section(~TF,:)]; %#ok<AGROW>
        end
        plot(reg(:,1),reg(:,2),'b.','LineWidth',1.5);
        plot(ol(:,1),ol(:,2),'r.','LineWidth',1.5);
        hold off
        return
    end
    
    fprintf('\nCleaning Up\n');
    fprintf('Allocated %d, Used %d\n', initLength, totalTrajUsed);
    % Remove extra space if available
    if totalTrajUsed < initLength 
        trajectories(totalTrajUsed+1:initLength) = [];
        trajPts(totalTrajUsed+1:initLength,:,:) = [];
        validPts(totalTrajUsed+1:initLength,:) = [];
    end
    
    
    % We're done with the frames, so mark live ones as done
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
%         if (trajectories{i}{5} == 1 || trajectories{i}{5} == 2) && trajectories{i}{3}(2)-trajectories{i}{3}(1) < p.Results.MinTrajLength
        if (trajectories{i}{5} == 1 || trajectories{i}{5} == 2) && trajectories{i}{3}(2)-trajectories{i}{3}(1) < 2
                trajectories{i}{5} = 0; % length 1 and 2 retired trajectories are invalid
                % Need at least length 3 to calculate acceleration for
                % optimization
                validPts(i,trajectories{i}{3}) = 0;
        elseif trajectories{i}{5} == 1
            trajectories{i}{5} = 2;
            trajectories{i}{3}(2) = ind-1;
        end
    end
    
    invalidInds = cellfun(@isempty,trajectories);
    trajPts(invalidInds,:,:) = [];
    validPts(invalidInds,:) = [];

    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        if trajectories{i}{5} == 0
            trajectories{i} = [];
            goodTraj = goodTraj - 1;
        end
    end
    
    S = whos('trajectories');
    fprintf('Valid Trajectories: %d\n', goodTraj);
    fprintf('mB: %f\n', S.bytes/1048576);
    
    if p.Results.ShowTrajectories
        figure
        hold on
        for i=1:numel(trajectories)
            if isempty(trajectories{i})
                continue
            end
            if trajectories{i}{5} == 0
                continue
            end
            pts = trajectories{i}{1};
            time = trajectories{i}{3};
            xpts = pts(:,1);
            ypts = pts(:,2);
            zpts = (time(1):time(2))';
            plot3(xpts, ypts, zpts, 'LineWidth', 1.5);
        end
        hold off
        view(0,0);
        return
    end
    
    %% Compute trajectory weights for each frame
    fprintf('\nCalculating Trajectory Weights\n');
    % numFrames x numTrajectories
    weights = zeros(numel(trajectories), numFrames);
    pd = makedist('Normal',0,p.Results.SigmaWeighting);
    numProcessed = 1;
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        if numProcessed ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*numProcessed/goodTraj, numProcessed, goodTraj);
        T = trajectories{i};
        ts = T{3}(1);
        te = T{3}(2);
        
        t = ts+1:te-1;
        ptsSeries = trajPts(numProcessed,:,:);
        diffs = trajPts(:,t,:)-ptsSeries(:,t,:);
        diffNorms = vecnorm(diffs,2,3);
        g = pdf(pd, diffNorms);
        g(validPts(:,t)<1) = 0;
        g(diffNorms>p.Results.SigmaWeighting*6) = 0;
        weights(i,t) = sum(g,1);
        ltmp = [t-ts; te-t];
        l2 = min(ltmp);
        l2(l2<0) = 0;
        weights(i,t) = l2 ./ weights(i,t);
        
        numProcessed = numProcessed + 1;
    end
    
    % Normalize weights
    fprintf('\n\nNormalizing\n');
    for i=1:numFrames
        if i ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*i/numFrames, i, numFrames);
        sumWeights = sum(weights(:,i));
        if sumWeights < 1e-12
            continue
        end
        weights(:,i) = weights(:,i) ./ sumWeights;
    end
    
    %% Optimization
    fprintf('\n\nOptimizing\n');
    % ld, ls, lu
    lambda = [1, 50, 0];
    f = @(x)VSObjectiveFun(x, numFrames, trajectories, weights, lambda, sz);
    % scale, theta, dx, dy
    initial = repmat([1.01; 0; 0; 0], numFrames, 1);
    lb = repmat([0.7; -Inf; -Inf; -Inf], numFrames, 1);
    ub = repmat([1.3; Inf; Inf; Inf], numFrames, 1);
%     options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
%         'SpecifyObjectiveGradient', true, 'ScaleProblem', 'jacobian', ...
%         'CheckGradients', false, 'FunctionTolerance', 1e-12, 'Display', 'iter');
%     options = optimoptions('fmincon', 'Algorithm', 'active-set', 'SpecifyObjectiveGradient', true, ...
%         'CheckGradients', false, 'MaxIterations', 1500, 'Display', 'iter');
    opt.Jacobian = 'on';
    opt.Display = 'iter';
    opt.DerivativeCheck = 'off';
    opt.MaxDamping = 1e12;
    opt.FactDamping = 5;
    tic
%     [T, resnorm, ~, exitflag, output] = lsqnonlin(f, initial, [], [], options);
%     [T, ~, exitflag, output] = fmincon(f, initial, [], [], [], [], lb, ub, [], options);
    [T, resnorm, ~, exitflag, ~, output] = LevenbergMarquardt(f, initial, lb, ub, opt);
    toc
    fprintf('Optimization Done\n');
%     fprintf('%s\n', output.message);
    fprintf('Squared Residual: %f\n', resnorm);
    fprintf('Exit Flag: %d\n', exitflag);
%     fprintf('Number of Iterations: %d\n', output.iterations);
%     fprintf('Number of function evaluations: %d\n', output.funcCount);
    
    stable = VideoWriter(output_name);
    open(stable);
    vid.CurrentTime = 0;
    index = 1;
    fprintf('\nStabilizing Video\n');
    while hasFrame(vid)
        if index ~= 1
            fprintf(repmat('\b', 1, nchar));
        end
        nchar = fprintf('%0.2f%% (%d/%d)', 100*index/numFrames, index, numFrames);
        frame = readFrame(vid);
        scale = T(4*(index-1)+1);
        theta = T(4*(index-1)+2);
        dx = T(4*(index-1)+3);
        dy = T(4*(index-1)+4);
        translation = [dx; dy];
        HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)], translation]; [0 0 1]];
        framewarp = imwarp(frame,affine2d(HsRt'),'OutputView',imref2d(size(frame)));
        writeVideo(stable, framewarp);
        index = index + 1;
        if index == numFrames+1
            break
        end
    end
    
    fprintf('\nDone\n');
    reshape(T, 4, [])'
    
    
    vid.CurrentTime = 0;
    frame = readFrame(vid);
    figure
    hold on
    imshow(frame);
    cnt = 1;
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        if trajectories{i}{3}(1) < 10
            pts = trajectories{i}{1};
            xpts = pts(:,1);
            ypts = pts(:,2);
            line(xpts, ypts, 'LineWidth', 1.5);
        end
        cnt = cnt + 1;
    end
    saveas(gcf, 'regular.png');
    
    figure
    hold on
    imshow(frame);
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        if trajectories{i}{3}(1) < 10
            pts = trajectories{i}{1};
            xptsT = zeros(size(pts,1),1);
            yptsT = xptsT;
            for j=1:size(pts,1)
                index = j+trajectories{i}{3}(1)-1;
                scale = T(4*(index-1)+1);
                theta = T(4*(index-1)+2);
                dx = T(4*(index-1)+3);
                dy = T(4*(index-1)+4);
                translation = [dx; dy];
                HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)], translation]; [0 0 1]];
                Tp = HsRt * [pts(j,:) 1]';
                xptsT(j) = Tp(1);
                yptsT(j) = Tp(2);
            end
            line(xptsT, yptsT, 'LineWidth', 1.5);
        end
        cnt = cnt + 1;
    end
    saveas(gcf, 'transformed.png');
end