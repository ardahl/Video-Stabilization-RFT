function StabilizeVideoRFTTest(varargin)
    sigmaSDefault = 15;
    p = inputParser;
    validPosNum = @(x) isnumeric(X) && isscalar(x) && x>=0;
    addOptional(p, 'sigmaWeighting', sigmaSDefault, validPosNum);
    parse(p, varargin{:});
    
    obj = TrajDebug();

    fprintf('\n%d Frames\n', obj.nFrames);
    
    % Set of all trajectories
    % Each trajectory is a cell of 5 elements
        % [x, y] - a rowx2 vector containing the x,y path
        % [features] - row vector containing the feature descriptions at
            % each time
        % [ts, te] - a 1x2 vector containing the start and end frames
        % [neighbor indices] - array of indices of the neighboring traj.
        % live = int on whether the trajectory is considered live(1),
            % retired(2), or invalid(0)

    
    % Process the first frame uniqely
    numFrames = obj.nFrames;
    fprintf('\nComputing Feature Trajectories\n');

    % Extract x,y from features
    pts = squeeze(obj.Points(1,:,:));
    % Allocate extra space so it doesn't resize as much
    % Assume average of length(features) features detected each frame
        % In order to reduce complexity at the expence of size, non-valid
        % trajectories of length 1 are kept.
    trajectories = cell(numFrames*size(pts,1), 1);
    trajPts = zeros(numFrames*size(pts,1), numFrames, 2);
    validPts = zeros(numFrames*size(pts,1), numFrames);
    totalTrajUsed = size(pts,1);
    initLength = numFrames*size(pts,1);
    % Compute Delaunay triangulation from points
    tri = delaunayTriangulation(pts);
    edge = edges(tri);
    % Record neighbor information for each feature
    neigh = cell(numFrames*size(pts,1),1);
    for i=1:size(edge,1)
        % Loop through each edge, need to mark down neighbor for (i,j) and
        % (j,i)
        e1 = edge(i,1);
        e2 = edge(i,2);
        neigh{e1} = [neigh{e1} e2];
        neigh{e2} = [neigh{e2} e1];
    end
    % Create new trajectories for each feature found
    for i=1:size(pts,1)
        newT = cell(1,5);
        newT{1} = pts(i,:);
        trajPts(i,1,:) = pts(i,:);
        validPts(i,1) = 1;
        newT{2} = 0;
        newT{3} = [1 1];
        newT{4} = neigh{i};
        newT{5} = 1;
        trajectories{i} = newT;
    end
    goodTraj = totalTrajUsed;
    fsize = 1;
    
    ind = 1;
    for frame=2:obj.nFrames

        %% Addition
        % Collect Features from the new frame
        
        %% Linking
        % Extract x,y from features
        pts = squeeze(obj.Points(frame,:,:));
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
        
        %% Propigation
        % Match features between the 2 images
        % Recover the features of all live trajectories
        trajFeat = zeros(goodTraj,fsize);
        trajInd = zeros(goodTraj,1);
        c1 = 1;
        for i=1:totalTrajUsed
            if isempty(trajectories{i})
                continue
            end
            if trajectories{i}{5} == 1
                trajFeat(c1) = trajectories{i}{2}(end,:);
                trajInd(c1) = i;
                c1 = c1 + 1;
            end
        end
        % IndexPairs(:,1) are the indices to the current trajectories 
        % (trajFeat) that are matching with the new features.
        % trajInd(indexPairs(:,1)) will give the indices to the trajectory
        % master list. They can then be updated with the indices to the new
        % features from indexPairs(:,2)
        indexPairs = obj.MatchInds(:,:,frame-1);
        pts2 = pts(indexPairs(:,2),:);
                
        % With the matching features, add them to the trajectories
        % Loop through each live trajectory and check if any of the
        % matching points from the previous from the previous frame are
        % part of the live set. 
        matchedTraj = trajInd(indexPairs(:,1));
        nodeInd = zeros(1,size(pts,1));
        for i=1:size(pts2,1)
            T = trajectories{matchedTraj(i)};
            newPts = [T{1}; pts2(i,:)];
            trajPts(matchedTraj(i),ind,:) = pts2(i,:);
            validPts(matchedTraj(i),ind) = 1;
            newFeat = [T{2}; 0];
            trajectories{matchedTraj(i)}{1} = newPts;
            trajectories{matchedTraj(i)}{2} = newFeat;
            trajectories{matchedTraj(i)}{3}(2) = ind;
            nodeInd(indexPairs(i,2)) = matchedTraj(i);
        end
        
        % For the new matches, create a new trajectory
        indMap = 1:size(pts,1);
        indMap(indexPairs(:,2)) = [];
        p2 = pts;
        p2(indexPairs(:,2),:) = [];
        for i=1:size(p2,1)
            newT = cell(1,5);
            newT{1} = p2(i,:);
            trajPts(totalTrajUsed+i,ind,:) = p2(i,:);
            validPts(totalTrajUsed+i,ind) = 1;
            newT{2} = 0;
            newT{3} = [ind ind];
            newT{4} = neigh{i};
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
        for i=1:size(pts,1)
            % Find trajectory that this feature belongs to, mark node i as
            % part of that trajectory
            if nodeInd(i) < 1
                fprintf('I dont think this should happen\n');
                continue
            end
            nList = neigh{i};
            nt = zeros(numel(nList));
            for j=1:numel(nList)
                nt(j) = nodeInd(j);
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
            if trajectories{i}{3}(2) ~= ind && trajectories{i}{3}(2)-trajectories{i}{3}(1) < 2
                trajectories{i}{5} = 0; % length 1 and 2 retired trajectories are invalid
                % Need at least length 3 to calculate acceleration for
                % optimization
                validPts(i,trajectories{i}{3}) = 0;
            elseif trajectories{i}{3}(2) ~= ind
                trajectories{i}{5} = 2;
            end
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
    end
    
    fprintf('\nCleaning Up\n');
    fprintf('Allocated %d, Used %d\n', initLength, totalTrajUsed);
    % Remove extra space if available
    if totalTrajUsed < initLength 
        trajectories(totalTrajUsed+1:initLength) = [];
%         neigh(totalTrajUsed+1:initLength) = [];
        trajPts(totalTrajUsed+1:initLength,:,:) = [];
        validPts(totalTrajUsed+1:initLength,:) = [];
    end
    
    % We're done with the frames, so mark live ones as done
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        if trajectories{i}{5} == 1 && trajectories{i}{3}(2)-trajectories{i}{3}(1) < 2
                trajectories{i}{5} = 0; % length 1 and 2 retired trajectories are invalid
                % Need at least length 3 to calculate acceleration for
                % optimization
                validPts(i,trajectories{i}{3}) = 0;
        elseif trajectories{i}{5} == 1
            trajectories{i}{5} = 2;
            trajectories{i}{3}(2) = ind;
        end
    end
    
    %s
    invalidInds = cellfun(@isempty,trajectories);
    trajPts(invalidInds,:,:) = [];
    validPts(invalidInds,:) = [];
    %e

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
    
    whos trajectories
    
    compareTraj(obj, trajectories);
    
    
    %% Compute trajectory weights for each frame
    fprintf('\nCalculating Trajectory Weights\n');
    % numFrames x numTrajectories
    weights = zeros(numel(trajectories), numFrames);
%     gaussNorm = 1/(p.Results.sigmaWeighting*sqrt(2*pi));
    pd = makedist('Normal',0,p.Results.sigmaWeighting);
    
    numProcessed = 1;
    for i=1:numel(trajectories)
        if isempty(trajectories{i})
            continue
        end
        T = trajectories{i};
        ts = T{3}(1);
        te = T{3}(2);
        
        t = ts+1:te-1;
        ptsSeries = trajPts(numProcessed,:,:);
        diffs = trajPts(:,t,:)-ptsSeries(:,t,:);
        diffNorms = vecnorm(diffs,2,3);
%         g = gaussNorm * exp((diffNorms.^2)./(2*(p.Results.sigmaWeighting).^2));
        g = pdf(pd, diffNorms);
        g(validPts(:,t)<1) = 0;
        g(diffNorms>p.Results.sigmaWeighting*6) = 0;
        weights(numProcessed,t) = sum(g,1);
        ltmp = [t-ts; te-t];
        l2 = min(ltmp);
        l2(l2<0) = 0;
        weights(i,t) = l2 ./ weights(i,t);
        
        numProcessed = numProcessed + 1;
    end

    % Normalize weights
    fprintf('\n\nNormalizing\n');
    for i=1:numFrames
        sumWeights = sum(weights(:,i));
        if sumWeights < 1e-12
            continue
        end
        weights(:,i) = weights(:,i) ./ sumWeights;
    end
    compareWeights(obj, weights);
    
    %% Optimization
    fprintf('\n\nOptimizing\n');
    % ld, ls, lu
    lambda = [1, 50, 1];
    sz = [800 800];
    f = @(x)VSObjectiveFun(x, numFrames, trajectories, weights, lambda, sz);
    Itrans = [1.01 0 0; 0 1.01 0; 0 0 1];
    initial = repmat(Itrans, numFrames, 1);
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
    [T, resnorm, ~, exitflag, output] = lsqnonlin(f, initial, [], [], options);
    fprintf('Optimization Done\n');
%     fprintf('%s\n', output.message);
    fprintf('Squared Residual: %f\n', resnorm);
    fprintf('Exit Flag: %d\n', exitflag);
    fprintf('Number of Iterations: %d\n', output.iterations);
    fprintf('Number of function evaluations: %d\n', output.funcCount);
    
    for i=1:numFrames
        T(i*3-2:i*3,:)
    end
    
    fprintf('Done\n');