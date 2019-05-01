function J = symJ(second)
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
Es = (s_-1)^4;
if exist('second','var') && second
    f = [Er, Es];
else
    f = Er;
end
J = jacobian(f, [sm_, thetam_, dxm_, dym_, s_, theta_, dx_, dy_, sp_, thetap_, dxp_, dyp_]);
% val = subs(J(1,:),[sm, thetam, dxm, dym, s, theta, dx, dy, sp, thetap, dxp, dyp, w, xm, ym, x, y, xp, yp], [1.01, 0, 0, 0, 1.01, 0, 0, 0, 1.01, 0, 0, 0, 1, 1, 1, 2, 2, 1, 1]);
simplify(J(4))
