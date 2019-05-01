s = 1;
t = 0;
dx = 0;
dy = 0;
scaleRand = 0.15;
% +/- 5 deg
% thetaRand = 0.0872665;
thetaRand = 0.26;
dxRand = 1;
dyRand = 1;

rng(1)
nFrames = 100;
k = 3;

Tparams = repmat([s, t, dx, dy], nFrames, 1);
base = Tparams + [scaleRand*(2*rand(nFrames,1)-1) thetaRand*(2*rand(nFrames,1)-1) dxRand*(2*rand(nFrames,1)-1) dyRand*(2*rand(nFrames,1)-1)];
baseT = zeros(3*nFrames,3);
smooth = base;
average = base;
% norms of (S-S2), (S-T), (S2-T)
% S - smoothed components matrix
% S2 - Averaged matrix
% T - Base
norms = zeros(nFrames,3);
for i=1:nFrames
    params = Tparams(i,:);
    s = params(1);
    t = params(2);
    tr = params(3:4);
    T = [[s*[cos(t) -sin(t); sin(t) cos(t)]; tr], [0 0 1]'];
    baseT(i*3-2:i*3,:) = T;
end

for i=1:nFrames
    % Base Matrix
    params = base(i,:);
    s = params(1);
    t = params(2);
    tr = params(3:4);
    T = [[s*[cos(t) -sin(t); sin(t) cos(t)]; tr], [0 0 1]'];
    
    % Average components
    lb = max(1, i-k);
    ub = min(nFrames, i+k);
    neighbors = base(lb:ub,:);
    ave = mean(neighbors);
    s = ave(1);
    t = ave(2);
    tr = ave(3:4);
    S = [[s*[cos(t) -sin(t); sin(t) cos(t)]; tr], [0 0 1]'];
    smooth(i,:) = [s t tr];
    
    S2 = zeros(3,3);
    for j=lb:ub
        s = base(j,1);
        t = base(j,2);
        tr = base(j,3:4);
        Stmp = [[s*[cos(t) -sin(t); sin(t) cos(t)]; tr], [0 0 1]'];
        S2 = S2 + Stmp;
    end
    S2 = S2 ./ (ub-lb+1);
    [s, t, tr] = decomposeT(S2);
    average(i,:) = [s t tr];
    
    norms(i,1) = norm(S-S2,'fro');
    norms(i,2) = norm(S-T,'fro');
    norms(i,3) = norm(S2-T,'fro');
end

fprintf('Max Error: %f\n', max(norms(:,1)));

% Plot everything
% Plot norms
figure
hold on
plot(norms(:,1), 'r-', 'LineWidth', 1.5);
% plot(norms(:,2), 'g-.', 'LineWidth', 1.5);
% plot(norms(:,3), 'b-.', 'LineWidth', 1.5);
title('Norms');
xlabel('Frame');
ylabel('Norm');
legend('S-S2','S-T','S2-T');
% Plot scale component
figure
subplot(1,2,1);
hold on
plot(base(:,1), 'g--', 'LineWidth', 1.5);
plot(smooth(:,1), 'r-', 'LineWidth', 1.5);
plot(average(:,1), 'b-', 'LineWidth', 1.5);
title('Scale');
xlabel('Frame');
ylabel('Scale');
legend('Base','Components','Matrices');
% Plot theta component
% figure
subplot(1,2,2);
hold on
plot(base(:,2), 'g--', 'LineWidth', 1.5);
plot(smooth(:,2), 'r-', 'LineWidth', 1.5);
plot(average(:,2), 'b-', 'LineWidth', 1.5);
title('Theta');
xlabel('Frame');
ylabel('Theta (rad)');
legend('Base','Components','Matrices');