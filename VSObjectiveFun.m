function [F, J] = VSObjectiveFun(x, nFrames, traj, weights, lambda, sz)
    % Form all the transformation matrices
    check_J = false;
    type = 1;
    % 1 for fi(x) to be each frame
    % 2 for fi(x) to be either Er or Edeg
    % 3 for 1 value for fmincon
    % 4 same as 1 but stretched by 4 (4*nFrames)
    if nargout > 1
        if type == 1 || type == 4
            J = zeros(nFrames,4*nFrames);
        else
            J = zeros(2, 4*nFrames);
        end
    end
    Tforms = zeros(3*nFrames,3);
    TformsInv = Tforms;
    for i=1:nFrames
        scale = x(4*(i-1)+1);
        scaleInv = 1/scale;
        theta = x(4*(i-1)+2);
        dx = x(4*(i-1)+3);
        dy = x(4*(i-1)+4);
        translation = [dx; dy];
        transInv = [-(dx*cos(theta)+dy*sin(theta)); (-dy*cos(theta)+dx*sin(theta))];
        Tforms(3*i-2:3*i,:) = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)], translation]; [0 0 1]];
        TformsInv(3*i-2:3*i,:) = [scaleInv*[[cos(theta) sin(theta); -sin(theta) cos(theta)], transInv]; [0 0 1]];
    end
    tol = 1e-12;
    F1 = zeros(nFrames,1);
    F2 = zeros(2,1);
    %% E_roughness
    Er = 0;
    for t=2:nFrames-1
        trajCount = 0;
        for i=1:numel(traj)
            if isempty(traj{i})
                continue
            end
            if t <= traj{i}{3}(1) || t >= traj{i}{3}(2)
                continue
            end
            w = weights(i, t);
            Tt = Tforms(t*3-2:t*3,:);
            Ttm1 = Tforms((t-1)*3-2:(t-1)*3,:);
            Ttp1 = Tforms((t+1)*3-2:(t+1)*3,:);
            pt = traj{i}{1}(t-traj{i}{3}(1)+1,:);
            ptm1 = traj{i}{1}((t-1)-traj{i}{3}(1)+1,:);
            ptp1 = traj{i}{1}((t+1)-traj{i}{3}(1)+1,:);
            a = Ttp1*[ptp1 1]' - 2*Tt*[pt 1]' + Ttm1*[ptm1 1]';
            Ttinv = TformsInv(t*3-2:t*3,:);
            normsq = norm(Ttinv*a).^2;
            val = w * normsq;
            Er = Er + val;
            F1(t) = F1(t) + val;
            if nargout > 1
                j = computeJEr(x, t, traj{i}, w, false);
                if check_J
                    j_check = computeJEr(x, t, traj{i}, w, true);
                    diff = j-j_check;
                    if max(diff,[],2) > 1e-6
                        fprintf('Difference Error %f\n', max(diff,2));
                        j
                        j_check
                        diff
                        return
                    end
                end
                if type == 1 || type == 4
                    J(t, 4*(t-2)+1:4*(t+1)) = J(t, 4*(t-2)+1:4*(t+1)) + j;
                else
                    J(1, 4*(t-2)+1:4*(t+1)) = J(1, 4*(t-2)+1:4*(t+1)) + j;
                end
            end
            trajCount = trajCount + 1;
        end
%         cutoff = 5;
%         if trajCount < cutoff
%             warning('Number of valid trajectories (%d) for frame %d is <%d. This may result in lower quality tranformations.',trajCount,t,cutoff);
%         end
    end
    %% E_degridation
    Edeg = 0;
    for t=1:nFrames
        %% E_distortion
        Ed = 0;
        
        % Turn the Transformation matrix into a similarity transform
        T = Tforms(t*3-2:t*3,:);
        Tinv = TformsInv(t*3-2:t*3,:);
        scale = x(4*(t-1)+1);
        %% E_scale
        Es = 0;
        if scale > 1
            Es = (scale-1).^4;
            if nargout > 1
                if type == 1 || type == 4
                    J(t,4*(t-1)+1) = J(t,4*(t-1)+1) + (4*(scale-1).^3);
                else
                    J(2,4*(t-1)+1) = 4*(scale-1).^3;
                end
            end
        end
        %% E_uncovered
        Eu = 0;
        delta = 0.001;
        v = [0 0 1; sz(1) 0 1; 0 sz(2) 1; sz(1) sz(2) 1];
        connEdges = [1 2; 2 4; 4 3; 3 1];
        for i=1:4
            v(i,:) = (Tinv*v(i,:)')';
            v(i,:) = v(i,:) / v(i,3);
        end
        e = [0 0 sz(1) 0; sz(1) 0 sz(1) sz(2); sz(1) sz(2) 0 sz(2); 0 sz(2) 0 0];
        center = [sz(1)/2 sz(2)/2];
        % Looping each edge
        for i=1:4
            edge = e(i,:);
            ref = (center(2)-edge(2))*(edge(3)-edge(1)) - (center(1)-edge(1))*(edge(4)-edge(2));
            len = norm(edge(3:4)-edge(1:2));
            lensq = len.^2;
            distDen = sqrt((edge(4)-edge(2)).^2 + (edge(3)-edge(1)).^2);
            % Looping each vertex
            for j=1:4
                if j ~= connEdges(i,1) && j ~= connEdges(i,2)
                    continue
                end
                vert = v(j,1:2);
                side = (vert(2)-edge(2))*(edge(3)-edge(1)) - (vert(1)-edge(1))*(edge(4)-edge(2));
                if sign(side) ~= 0 && sign(ref) == sign(side)
                    pe1 = vert - edge(1:2);
                    e1e2 = edge(3:4) - edge(1:2);
                    dot = pe1(1)*e1e2(1)+pe1(2)*e1e2(2);
                    time = max(0, min(1, dot/lensq));
                    proj = edge(1:2) + time * (edge(3:4) - edge(1:2));
                    dist = vert - proj;
                    dist = norm(dist);
                    
                    dist = sqrt(dist.^2 + delta.^2) - delta;
                    Eu = Eu + (len*dist).^2;
                end
            end
        end
        
        Edeg = Edeg + (Ed + lambda(2)*Es + lambda(3)*Eu);
        F1(t) = F1(t) + lambda(1)*(Ed + lambda(2)*Es + lambda(3)*Eu);
    end
    F2(1) = Er;
    F2(2) = lambda(1)*Edeg;
    
    if type == 1
        F = F1;
    elseif type == 2
        F = F2;
    elseif type == 3
        F = F2(1) + F2(2);
        if nargout > 1
            J = sum(J,1);
        end
    elseif type == 4
        F = repelem(F1, 4);
        J = repelem(J, 4, 1);
    end
end


function J = computeJEr(X, t, traj, w, useSym)
    sm = X(4*(t-2)+1);
    thetam = X(4*(t-2)+2);
    dxm = X(4*(t-2)+3);
    dym = X(4*(t-2)+4);
    ptm = traj{1}((t-1)-traj{3}(1)+1,:);
    xm = ptm(1);
    ym = ptm(2);
    
    s = X(4*(t-1)+1);
    theta = X(4*(t-1)+2);
    dx = X(4*(t-1)+3);
    dy = X(4*(t-1)+4);
    pt = traj{1}(t-traj{3}(1)+1,:);
    x = pt(1);
    y = pt(2);
    
    sp = X(4*t+1);
    thetap = X(4*t+2);
    dxp = X(4*t+3);
    dyp = X(4*t+4);
    ptp = traj{1}((t+1)-traj{3}(1)+1,:);
    xp = ptp(1);
    yp = ptp(2);
    
    if useSym
        syms s_ theta_ dx_ dy_
        syms sm_ thetam_ dxm_ dym_
        syms sp_ thetap_ dxp_ dyp_
        syms x_ y_ xm_ ym_ xp_ yp_
        syms w_
        T = [s_*cos(theta_) -s_*sin(theta_) dx_;
             s_*sin(theta_)  s_*cos(theta_) dy_;
                   0               0          1];
        Tm = [sm_*cos(thetam_) -sm_*sin(thetam_) dxm_;
              sm_*sin(thetam_)  sm_*cos(thetam_) dym_;
                     0                 0            1];
        Tp = [sp_*cos(thetap_) -sp_*sin(thetap_) dxp_;
              sp_*sin(thetap_)  sp_*cos(thetap_) dyp_;
                     0                 0            1];
        Tinv = [cos(theta_)/s_  sin(theta_)/s_ -(dx_*cos(theta_)+dy_*sin(theta_))/s_;
                -sin(theta_)/s_ cos(theta_)/s_ (-dy_*cos(theta_)+dx_*sin(theta_))/s_;
                      0                0                         1                  ];
        Pt = [x_ y_ 1]';
        Pm = [xm_ ym_ 1]';
        Pp = [xp_ yp_ 1]';
        a = Tp*Pp - 2*T*Pt + Tm*Pm;
        nrm = Tinv*a;
        Er = w_ * (nrm(1)^2 + nrm(2)^2);
        Jac = jacobian(Er, [sm_, thetam_, dxm_, dym_, s_, theta_, dx_, dy_, sp_, thetap_, dxp_, dyp_]);
        val = subs(Jac,[sm_, thetam_, dxm_, dym_, s_, theta_, dx_, dy_, sp_, thetap_, dxp_, dyp_, w_, xm_, ym_, x_, y_, xp_, yp_], ...
            [sm, thetam, dxm, dym, s, theta, dx, dy, sp, thetap, dxp, dyp, w, xm, ym, x, y, xp, yp]);
        J = double(val);
        return
    end
    
    sinThetam = sin(thetam);
    cosThetam = cos(thetam);
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    sinThetap = sin(thetap);
    cosThetap = cos(thetap);
    
    cosThetatm = cos(theta-thetam);
    cosThetamp = cos(thetam-thetap);
    cosThetatp = cos(theta-thetap);
    sinThetatm = sin(theta-thetam);
    sinThetamp = sin(thetam-thetap);
    sinThetatp = sin(theta-thetap);
    
    
    Dsm = 2 * (-ym*cosThetatm+xm*sinThetatm) * (2*s*y-(-2*dy+dym+dyp)*cosTheta ...
          - sm*ym*cosThetatm - sp*yp*cosThetatp - 2*dx*sinTheta + dxm*sinTheta ...
          + dxp*sinTheta + sm*xm*sinThetatm + sp*xp*sinThetatp) + 2*(xm*cosThetatm+ym*sinThetatm) ...
          * (-2*s*x + (-2*dx+dxm+dxp)*cosTheta + sm*xm*cosThetatm + sp*xp*cosThetatp ...
          - 2*dy*sinTheta + dym*sinTheta + dyp*sinTheta + sm*ym*sinThetatm + sp*yp*sinThetatp);
    Dsm = Dsm * (w/s.^2);
    Dtm = 2*s*(xm*y-x*ym)*cosThetatm + (2*dy*xm-dym*xm-dyp*xm-2*dx*ym+dxm*ym+dxp*ym)*cosThetam ...
          + sp*xp*ym*cosThetamp - sp*xm*yp*cosThetamp + 2*s*x*xm*sinThetatm ...
          + 2*s*y*ym*sinThetatm - 2*dx*xm*sinThetam + dxm*xm*sinThetam ...
          + dxp*xm*sinThetam - 2*dy*ym*sinThetam + dym*ym*sinThetam ...
          + dyp*ym*sinThetam + sp*xm*xp*sinThetamp + sp*ym*yp*sinThetamp;
    Dtm = Dtm * (-(2*w*sm)/s.^2);
    Ddxm = -2*dx + dxm + dxp - 2*s*x*cosTheta + sm*xm*cosThetam ...
           + sp*xp*cosThetap + 2*s*y*sinTheta - sm*ym*sinThetam - sp*yp*sinThetap;
    Ddxm = Ddxm * ((2*w)/s.^2);
    Ddym = -2*dy + dym + dyp - 2*s*y*cosTheta + sm*ym*cosThetam ...
           + sp*yp*cosThetap - 2*s*x*sinTheta + sm*xm*sinThetam + sp*xp*sinThetap;
    Ddym = Ddym * ((2*w)/s.^2);
    
    
    Ds = 4*s*y*(2*s*y-(-2*dy+dym+dyp)*cosTheta-sm*ym*cosThetatm-sp*yp*cosThetatp-2*dx*sinTheta+dxm*sinTheta+dxp*sinTheta+sm*xm*sinThetatm+sp*xp*sinThetatp) ...
         - 4*s*x*(-2*s*x+(-2*dx+dxm+dxp)*cosTheta+sm*xm*cosThetatm+sp*xp*cosThetatp-2*dy*sinTheta+dym*sinTheta+dyp*sinTheta+sm*ym*sinThetatm+sp*yp*sinThetatp) ...
         - 2*((2*s*y-(-2*dy+dym+dyp)*cosTheta-sm*ym*cosThetatm-sp*yp*cosThetatp-2*dx*sinTheta+dxm*sinTheta+dxp*sinTheta+sm*xm*sinThetatm+sp*xp*sinThetatp).^2 ...
         + (-2*s*x+(-2*dx+dxm+dxp)*cosTheta+sm*xm*cosThetatm+sp*xp*cosThetatp-2*dy*sinTheta+dym*sinTheta+dyp*sinTheta+sm*ym*sinThetatm+sp*yp*sinThetatp).^2);
    Ds = Ds * (w/s.^3);
    Dt = (2*dy*x-dym*x-dyp*x-2*dx*y+dxm*y+dxp*y)*cosTheta + sm*(xm*y-x*ym)*cosThetatm ...
         + sp*xp*y*cosThetatp - sp*x*yp*cosThetatp - 2*dx*x*sinTheta + dxm*x*sinTheta ...
         + dxp*x*sinTheta - 2*dy*y*sinTheta + dym*y*sinTheta + dyp*y*sinTheta ...
         + sm*x*xm*sinThetatm + sm*y*ym*sinThetatm + sp*x*xp*sinThetatp + sp*y*yp*sinThetatp;
    Dt = Dt * ((4*w)/s);
    Ddx = 2*dx - dxm - dxp + 2*s*x*cosTheta - sm*xm*cosThetam - sp*xp*cosThetap ...
          - 2*s*y*sinTheta + sm*ym*sinThetam + sp*yp*sinThetap;
    Ddx = Ddx * ((4*w)/s.^2);
    Ddy = 2*dy - dym - dyp + 2*s*y*cosTheta - sm*ym*cosThetam - sp*yp*cosThetap ...
          + 2*s*x*sinTheta - sm*xm*sinThetam - sp*xp*sinThetap;
    Ddy = Ddy * ((4*w)/s.^2);
    
    
    Dsp = 2 * (-yp*cosThetatp+xp*sinThetatp) * (2*s*y-(-2*dy+dym+dyp)*cosTheta ...
          - sm*ym*cosThetatm - sp*yp*cosThetatp - 2*dx*sinTheta + dxm*sinTheta ...
          + dxp*sinTheta + sm*xm*sinThetatm + sp*xp*sinThetatp) + 2*(xp*cosThetatp+yp*sinThetatp) ...
          * (-2*s*x + (-2*dx+dxm+dxp)*cosTheta + sm*xm*cosThetatm + sp*xp*cosThetatp ...
          - 2*dy*sinTheta + dym*sinTheta + dyp*sinTheta + sm*ym*sinThetatm + sp*yp*sinThetatp);
    Dsp = Dsp * (w/s.^2);
    Dtp = 2*s*(xp*y-x*yp)*cosThetatp + sm*(-xp*ym+xm*yp)*cosThetamp + 2*dy*xp*cosThetap ...
          - dym*xp*cosThetap - dyp*xp*cosThetap - 2*dx*yp*cosThetap + dxm*yp*cosThetap ...
          + dxp*yp*cosThetap + 2*s*x*xp*sinThetatp + 2*s*y*yp*sinThetatp - sm*xm*xp*sinThetamp ...
          - sm*ym*yp*sinThetamp - 2*dx*xp*sinThetap + dxm*xp*sinThetap + dxp*xp*sinThetap ...
          - 2*dy*yp*sinThetap + dym*yp*sinThetap + dyp*yp*sinThetap;
    Dtp = Dtp * (-(2*w*sp)/s.^2);
    Ddxp = Ddxm;
    Ddyp = Ddym;
    
    J = [Dsm Dtm Ddxm Ddym Ds Dt Ddx Ddy Dsp Dtp Ddxp Ddyp];
end