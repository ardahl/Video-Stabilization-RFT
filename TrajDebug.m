classdef TrajDebug
    % Contains enough information to use for debugging
    %   List of points of features for each frame (nFramesxnPtsx2)
    %   List of matching indices ind where the ind(feat,i)th point in frame
    %       t matches with the ind(feat,j)th point in t+1
    %
    properties
        nFrames=5;
        nPts=7;
        sigma=15;
        Points
        MatchInds
        ActualTrajs
        ActualWeights
    end
    methods
        function obj=TrajDebug()
            Points = zeros(obj.nFrames, obj.nPts, 2);
            diff = 20;
%             Points(1,1,:) = [1 1];
%             Points(1,2,:) = [500 132];
%             Points(1,3,:) = [416 210];
%             Points(1,4,:) = [91 384];
%             Points(1,5,:) = [250 250];
            Points(1,1,:) = [200 200];
            Points(1,2,:) = [250 250];
            Points(1,3,:) = [300 300];
            Points(1,4,:) = [200 300];
            Points(1,5,:) = [300 200];
            Points(1,6,:) = [200 324];
            Points(1,7,:) = [37 98];
            Points(2,1,:) = diff+Points(1,1,:);
            Points(2,2,:) = diff+Points(1,2,:);
            Points(2,3,:) = diff+Points(1,3,:);
            Points(2,4,:) = diff+Points(1,4,:);
            Points(2,5,:) = diff+Points(1,5,:);
            Points(2,6,:) = [17 192];
            Points(2,7,:) = [198 412];
            Points(3,1,:) = diff+Points(2,1,:);
            Points(3,2,:) = diff+Points(2,2,:);
            Points(3,3,:) = diff+Points(2,3,:);
            Points(3,4,:) = diff+Points(2,4,:);
            Points(3,5,:) = diff+Points(2,5,:);
            Points(3,6,:) = [277 391];
            Points(3,7,:) = [200 210];
            Points(4,1,:) = -diff+Points(3,1,:);
            Points(4,2,:) = -diff+Points(3,2,:);
            Points(4,3,:) = -diff+Points(3,3,:);
            Points(4,4,:) = -diff+Points(3,4,:);
            Points(4,5,:) = -diff+Points(3,5,:);
            Points(4,6,:) = [10 10];
            Points(4,7,:) = [11 11];
            Points(5,1,:) = -diff+Points(4,1,:);
            Points(5,2,:) = -diff+Points(4,2,:);
            Points(5,3,:) = -diff+Points(4,3,:);
            Points(5,4,:) = -diff+Points(4,4,:);
            Points(5,5,:) = -diff+Points(4,5,:);
            Points(5,6,:) = [12 12];
            Points(5,7,:) = [13 13];
            
            MatchInds = zeros(5,2,obj.nFrames-1);
            MatchInds(1,:,1) = [1 1];
            MatchInds(2,:,1) = [2 2];
            MatchInds(3,:,1) = [3 3];
            MatchInds(4,:,1) = [4 4];
            MatchInds(5,:,1) = [5 5];
            MatchInds(1,:,2) = [1 1];
            MatchInds(2,:,2) = [2 2];
            MatchInds(3,:,2) = [3 3];
            MatchInds(4,:,2) = [4 4];
            MatchInds(5,:,2) = [5 5];
            MatchInds(1,:,3) = [1 1];
            MatchInds(2,:,3) = [2 2];
            MatchInds(3,:,3) = [3 3];
            MatchInds(4,:,3) = [4 4];
            MatchInds(5,:,3) = [5 5];
            MatchInds(1,:,4) = [1 1];
            MatchInds(2,:,4) = [2 2];
            MatchInds(3,:,4) = [3 3];
            MatchInds(4,:,4) = [4 4];
            MatchInds(5,:,4) = [5 5];
            
            obj.Points = Points;
            obj.MatchInds = MatchInds;
            
            % Fill out what the trajectories should be
            ActualTraj = cell(5,1);
            % First frame
            newT = cell(1,5);
            newT{1} = Points(1,1,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{1} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,2,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{2} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,3,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{3} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,4,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{4} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,5,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{5} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,6,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{6} = newT;
            newT = cell(1,5);
            newT{1} = Points(1,7,:);
            newT{2} = 0;
            newT{3} = [1 1];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{7} = newT;
            % Second Frame
            T = ActualTraj{1};
            newPts = [T{1}; Points(2,1,:)];
            newFeat = [T{2}; 0];
            ActualTraj{1}{1} = newPts;
            ActualTraj{1}{2} = newFeat;
            ActualTraj{1}{3}(2) = 2;
            T = ActualTraj{2};
            newPts = [T{1}; Points(2,2,:)];
            newFeat = [T{2}; 0];
            ActualTraj{2}{1} = newPts;
            ActualTraj{2}{2} = newFeat;
            ActualTraj{2}{3}(2) = 2;
            T = ActualTraj{3};
            newPts = [T{1}; Points(2,3,:)];
            newFeat = [T{2}; 0];
            ActualTraj{3}{1} = newPts;
            ActualTraj{3}{2} = newFeat;
            ActualTraj{3}{3}(2) = 2;
            T = ActualTraj{4};
            newPts = [T{1}; Points(2,4,:)];
            newFeat = [T{2}; 0];
            ActualTraj{4}{1} = newPts;
            ActualTraj{4}{2} = newFeat;
            ActualTraj{4}{3}(2) = 2;
            T = ActualTraj{5};
            newPts = [T{1}; Points(2,5,:)];
            newFeat = [T{2}; 0];
            ActualTraj{5}{1} = newPts;
            ActualTraj{5}{2} = newFeat;
            ActualTraj{5}{3}(2) = 2;
            newT = cell(1,5);
            newT{1} = Points(2,6,:);
            newT{2} = 0;
            newT{3} = [2 2];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{8} = newT;
            newT = cell(1,5);
            newT{1} = Points(2,7,:);
            newT{2} = 0;
            newT{3} = [2 2];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{9} = newT;
            % Third Frame
            T = ActualTraj{1};
            newPts = [T{1}; Points(3,1,:)];
            newFeat = [T{2}; 0];
            ActualTraj{1}{1} = newPts;
            ActualTraj{1}{2} = newFeat;
            ActualTraj{1}{3}(2) = 3;
            T = ActualTraj{2};
            newPts = [T{1}; Points(3,2,:)];
            newFeat = [T{2}; 0];
            ActualTraj{2}{1} = newPts;
            ActualTraj{2}{2} = newFeat;
            ActualTraj{2}{3}(2) = 3;
            T = ActualTraj{3};
            newPts = [T{1}; Points(3,3,:)];
            newFeat = [T{2}; 0];
            ActualTraj{3}{1} = newPts;
            ActualTraj{3}{2} = newFeat;
            ActualTraj{3}{3}(2) = 3;
            T = ActualTraj{4};
            newPts = [T{1}; Points(3,4,:)];
            newFeat = [T{2}; 0];
            ActualTraj{4}{1} = newPts;
            ActualTraj{4}{2} = newFeat;
            ActualTraj{4}{3}(2) = 3;
            T = ActualTraj{5};
            newPts = [T{1}; Points(3,5,:)];
            newFeat = [T{2}; 0];
            ActualTraj{5}{1} = newPts;
            ActualTraj{5}{2} = newFeat;
            ActualTraj{5}{3}(2) = 3;
            newT = cell(1,5);
            newT{1} = Points(3,6,:);
            newT{2} = 0;
            newT{3} = [3 3];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{10} = newT;
            newT = cell(1,5);
            newT{1} = Points(3,7,:);
            newT{2} = 0;
            newT{3} = [3 3];
            newT{4} = 0;
            newT{5} = 1;
            ActualTraj{11} = newT;
            % Fourth Frame
            T = ActualTraj{1};
            newPts = [T{1}; Points(4,1,:)];
            newFeat = [T{2}; 0];
            ActualTraj{1}{1} = newPts;
            ActualTraj{1}{2} = newFeat;
            ActualTraj{1}{3}(2) = 4;
            T = ActualTraj{2};
            newPts = [T{1}; Points(4,2,:)];
            newFeat = [T{2}; 0];
            ActualTraj{2}{1} = newPts;
            ActualTraj{2}{2} = newFeat;
            ActualTraj{2}{3}(2) = 4;
            T = ActualTraj{3};
            newPts = [T{1}; Points(4,3,:)];
            newFeat = [T{2}; 0];
            ActualTraj{3}{1} = newPts;
            ActualTraj{3}{2} = newFeat;
            ActualTraj{3}{3}(2) = 4;
            T = ActualTraj{4};
            newPts = [T{1}; Points(4,4,:)];
            newFeat = [T{2}; 0];
            ActualTraj{4}{1} = newPts;
            ActualTraj{4}{2} = newFeat;
            ActualTraj{4}{3}(2) = 4;
            T = ActualTraj{5};
            newPts = [T{1}; Points(4,5,:)];
            newFeat = [T{2}; 0];
            ActualTraj{5}{1} = newPts;
            ActualTraj{5}{2} = newFeat;
            ActualTraj{5}{3}(2) = 4;
            newT = cell(1,5);
            newT{1} = Points(4,6,:);
            newT{2} = 0;
            newT{3} = [4 4];
            newT{4} = 0;
            newT{5} = 0;
            ActualTraj{12} = newT;
            newT = cell(1,5);
            newT{1} = Points(4,7,:);
            newT{2} = 0;
            newT{3} = [4 4];
            newT{4} = 0;
            newT{5} = 0;
            ActualTraj{13} = newT;
            % Fifth Frame
            T = ActualTraj{1};
            newPts = [T{1}; Points(5,1,:)];
            newFeat = [T{2}; 0];
            ActualTraj{1}{1} = newPts;
            ActualTraj{1}{2} = newFeat;
            ActualTraj{1}{3}(2) = 5;
            ActualTraj{1}{5} = 2;
            T = ActualTraj{2};
            newPts = [T{1}; Points(5,2,:)];
            newFeat = [T{2}; 0];
            ActualTraj{2}{1} = newPts;
            ActualTraj{2}{2} = newFeat;
            ActualTraj{2}{3}(2) = 5;
            ActualTraj{2}{5} = 2;
            T = ActualTraj{3};
            newPts = [T{1}; Points(5,3,:)];
            newFeat = [T{2}; 0];
            ActualTraj{3}{1} = newPts;
            ActualTraj{3}{2} = newFeat;
            ActualTraj{3}{3}(2) = 5;
            ActualTraj{3}{5} = 2;
            T = ActualTraj{4};
            newPts = [T{1}; Points(5,4,:)];
            newFeat = [T{2}; 0];
            ActualTraj{4}{1} = newPts;
            ActualTraj{4}{2} = newFeat;
            ActualTraj{4}{3}(2) = 5;
            ActualTraj{4}{5} = 2;
            T = ActualTraj{5};
            newPts = [T{1}; Points(5,5,:)];
            newFeat = [T{2}; 0];
            ActualTraj{5}{1} = newPts;
            ActualTraj{5}{2} = newFeat;
            ActualTraj{5}{3}(2) = 5;
            ActualTraj{5}{5} = 2;
            newT = cell(1,5);
            newT{1} = Points(5,6,:);
            newT{2} = 0;
            newT{3} = [5 5];
            newT{4} = 0;
            newT{5} = 0;
            ActualTraj{14} = newT;
            newT = cell(1,5);
            newT{1} = Points(5,7,:);
            newT{2} = 0;
            newT{3} = [5 5];
            newT{4} = 0;
            newT{5} = 0;
            ActualTraj{15} = newT;
            
            for i=6:15
                ActualTraj{i} = [];
            end
            obj.ActualTrajs = ActualTraj;
            
            obj.ActualWeights = zeros(15, 5);
            diffNorm = zeros(15,5);
            pd = makedist('Normal',0,obj.sigma);
            for i=1:5 %Traj
                for j=1:5 %Frames
                    obj.ActualWeights(i,j) = min(j-1,5-j);
                    Pi = ActualTraj{i}{1}(j,:);
                    for k=1:5
                        Pj = ActualTraj{i}{1}(k,:);
                        dn = norm(Pj-Pi);
                        if dn < obj.sigma*6
                            diffNorm(i,j) = diffNorm(i,j) + pdf(pd,dn);
                        end
                    end
                end
            end
            obj.ActualWeights
            diffNorm
            for i=1:5 %Traj
                for j=1:5 %frames
                    obj.ActualWeights(i,j) = obj.ActualWeights(i,j) / diffNorm(i,j);
                end
            end
            sumWeights = sum(obj.ActualWeights, 1);
            obj.ActualWeights
            sumWeights
            for i=1:5 %Traj
                for j=1:5 %frames
                    if sumWeights(j) > 1e-12
                        obj.ActualWeights(i,j) = obj.ActualWeights(i,j) / sumWeights(j);
                    end
                end
            end
        end
        function compareTraj(obj, Traj)
            fprintf('Sizes: %d, %d\n', numel(obj.ActualTrajs), numel(Traj));
            if numel(obj.ActualTrajs) ~= numel(Traj)
                fprintf('Different Base Sizes\n');
                return
            end
            valid1 = 0;
            valid2 = 0;
            for i=1:numel(obj.ActualTrajs)
                if ~isempty(obj.ActualTrajs{i})
                    valid1 = valid1 + 1;
                end
            end
            for i=1:numel(obj.ActualTrajs)
                if ~isempty(Traj{i})
                    valid2 = valid2 + 1;
                end
            end
            fprintf('Valid: %d, %d\n', valid1, valid2);
            if valid1 ~= valid2
                fprintf('Different Valid Sizes\n');
                return
            end
            for i=1:numel(obj.ActualTrajs)
                if (isempty(obj.ActualTrajs{i}) && ~isempty(Traj{i})) || (~isempty(obj.ActualTrajs{i}) && isempty(Traj{i}))
                    fprintf('Trajectory %d Different Order\n', i);
                    return
                end
            end
            for i=1:numel(obj.ActualTrajs)
                if isempty(obj.ActualTrajs{i})
                    continue
                end
                % Compare Points
                p1 = obj.ActualTrajs{i}{1};
                p2 = Traj{i}{1};
                if p1 ~= p2
                    p1 %#ok<*NOPRT>
                    p2
                    fprintf('Trajectory %d Different Path\n', i);
                    return
                end
                % Compare timestamps
                ts1 = obj.ActualTrajs{i}{3};
                ts2 = Traj{i}{3};
                if ts1 ~= ts2
                    ts1
                    ts2
                    fprintf('Trajectory %d Different Timestamp\n', i);
                    return
                end
                % Compare flag
                f1 = obj.ActualTrajs{i}{5};
                f2 = Traj{i}{5};
                if f1 ~= f2
                    f1
                    f2
                    fprintf('Trajectory %d Different Flags\n', i);
                    return
                end
            end
            fprintf('Trajectories are equal\n');
        end
        function compareWeights(obj, weights)
            fprintf('Sizes: (%d,%d), (%d,%d)\n', size(obj.ActualWeights,1), size(obj.ActualWeights,2), size(weights,1), size(weights,2));
            if numel(obj.ActualWeights) ~= numel(weights)
                fprintf('Different Base Sizes\n');
                return
            end
            if any(any(abs(obj.ActualWeights-weights)>1e-3))
                obj.ActualWeights
                weights
                abs(obj.ActualWeights-weights) > 1e-12
                fprintf('Different Weighting Matrices\n');
                return
            end
            fprintf('Weights Equal\n');
        end
    end
end