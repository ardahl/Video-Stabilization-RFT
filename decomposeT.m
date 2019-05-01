function [s, theta, translation] = decomposeT(T)
    R = T(1:2,1:2);
    % Compute theta from mean of two possible arctangents
    theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
    % Compute scale from mean of two stable mean calculations
    s = mean(R([1 4])/cos(theta));
    % Translation remains the same:
    translation = T(3, 1:2);
end