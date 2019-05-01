function [c, ceq] = scaleConstraint(x, numFrames)
    c = zeros(numFrames,1);
    for t=1:numFrames
        T = x(t*3-2:t*3,:);
        R = T(1:2,1:2);
        % Compute theta from mean of two possible arctangents
        theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
        % Compute scale from mean of two stable mean calculations
        scale = mean(R([1 4])/cos(theta));
        c(t) = -scale;
    end
    ceq = [];
end