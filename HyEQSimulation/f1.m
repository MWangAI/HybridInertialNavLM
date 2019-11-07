function xdot = ffun(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  A1 A2 g e3 sigmaa sigmag

%% state
% xIn = [x0;xhat0;vP10;xhat0;vP20;tau;time];  % 15+(15+3)+(15+45) +2
% real system states
vR     = x(1:9,1);
p      = x(10:12,1);
v      = x(13:15,1);
% A-HINO states 
index  = 15;
lP1    = length(A1)*(length(A1)+1)/2;
vRhat1 = x(index+1:index+9,1);
phat1  = x(index+10:index+12,1);
vhat1  = x(index+13:index+15,1);
vP1    = x(index+16:index+15+lP1,1);
% iEKF states
index  = index+15+lP1;
lP2    = length(A2)*(length(A2)+1)/2;
vRhat2 = x(index+1:index+9,1);
phat2  = x(index+10:index+12,1);
vhat2  = x(index+13:index+15,1);
vP2    = x(index+16:index+15+lP2,1);
% time
tau    = x(end-1,1);
time   = x(end,1); 
 
% omega and accel measurements
R         = reshape(vR,3,3);
[a,omega] = imu(R,time);
a_y       = a + sigmaa*randn(3,1);
omega_y   = omega + sigmag*randn(3,1);

% dynamics of real system 
DR    = R*Skew(omega);
Dp    = v;
Dv    = -g*e3 + R*a;
DvR   = reshape(DR,9,1);
Dx    = [DvR;Dp;Dv];

% A-HINO 
Rhat1   = reshape(vRhat1,3,3);
DRhat1  = Rhat1*Skew(omega_y);
Dphat1  = vhat1;
Dvhat1  = -g*e3 + Rhat1*a_y;
DvRhat1 = reshape(DRhat1,9,1);
P1      = reshapeT(vP1);
V1      = 0.1*mean(sigmaa.^2)*diag([0.001 1]);
% V1      = diag([0.0*ones(1,3) sigmaa.^2]);
DP1     = A1*P1+P1*A1' + V1;
DvP1    = reshapeT(DP1);
Dxhat1  = [DvRhat1;Dphat1;Dvhat1;DvP1];

% iEKF  
Rhat2   = reshape(vRhat2,3,3);
Xhat2   = [Rhat2 vhat2 phat2;zeros(2,3) eye(2)];
P2      = reshapeT(vP2);
DRhat2  = Rhat2*Skew(omega_y);
Dphat2  = vhat2;
Dvhat2  = -g*e3 + Rhat2*a_y;
DvRhat2 = reshape(DRhat2,9,1);
P2      = reshapeT(vP2);
CovW    = diag([sigmag.^2, sigmaa.^2 zeros(1,3)]);
Temp    = [Rhat2 zeros(3,6);Skew(vhat2)*Rhat2 Rhat2 zeros(3);
           Skew(phat2)*Rhat2 zeros(3) Rhat2];
V2      = Temp*CovW*Temp';
DP2     = A2*P2+P2*A2' + V2;
DvP2    = reshapeT(DP2);
Dxhat2  = [DvRhat2;Dphat2;Dvhat2;DvP2];


% output
xdot=[Dx;Dxhat1;Dxhat2;-1;1];

end