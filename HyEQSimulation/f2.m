function xplus = gfun(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A1 A2 C1 C2 pointsNum   Kn pointsI  sigmap kR

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
 
%% omega and accel measurements
R       = reshape(vR,3,3);
pointsB = R'*(pointsI - p) + sigmap*randn(size(pointsI));
r       = [pointsI;zeros(1,pointsNum);ones(1,pointsNum)];
b       = [pointsB;zeros(1,pointsNum);ones(1,pointsNum)];

% dynamics of real system 
vx      = [vR;p;v];

% A-HINO 
Rhat1   = reshape(vRhat1,3,3);
Xhat1   = [Rhat1 vhat1 phat1;zeros(2,3) eye(2)];
pointsC = sum(pointsI*Kn,2)/pointsNum; % center of landmarks 
Xc      = [eye(3) zeros(3,1) pointsC;zeros(2,3) eye(2)];
XcInv   = [eye(3) zeros(3,1) -pointsC;zeros(2,3) eye(2)];
Xhat1   = RestX(Xhat1,r,b,Kn);
P1      = reshapeT(vP1);
Q1      = 1*sigmap^2;    %1.0e-4
Kvp     = P1*C1'/(C1*P1*C1'+Q1);
P1      = (eye(size(P1))-Kvp*C1)*P1;
kc      = sum(diag(Kn));
% Temp    = projectionSE2(XcInv*(r-Xhat1*b)*Kn*r'*XcInv',kR,Kvp);
K       = [kR*eye(3) zeros(3,2);zeros(1,5);zeros(1,3) Kvp(2,1)/kc Kvp(1,1)/kc];
Temp    = projectionSE(XcInv*(r-Xhat1*b)*Kn*r'*XcInv'*K);
Delta1  = -Xc*Temp*XcInv;
Xhat1   = expm(-Delta1)*Xhat1;

vRhat1  = reshape(Xhat1(1:3,1:3),9,1);
vhat1   = Xhat1(1:3,4);
phat1   = Xhat1(1:3,5);
vP1     = reshapeT(P1);
vxhat1  = [vRhat1;phat1;vhat1;vP1];

% iEKF  
Rhat2   = reshape(vRhat2,3,3);
Xhat2   = [Rhat2 vhat2 phat2;zeros(2,3) eye(2)];
P2      = reshapeT(vP2);
Q2      = sigmap^2*eye(3*pointsNum);
Ln      = P2*C2'/(C2*P2*C2'+Q2);
P2      = (eye(size(P2))-Ln*C2)*P2;
II      = [eye(3) zeros(3,2)];
Temp    = II*(Xhat2*b-r);
Temp    = Ln*reshape(Temp,3*pointsNum,1);
Delta2  = [Skew(Temp(1:3,1)) Temp(4:6,1) Temp(7:9,1);zeros(2,5)];
Xhat2   = expm(Delta2)*Xhat2;

vRhat2  = reshape(Xhat2(1:3,1:3),9,1);
vhat2   = Xhat2(1:3,4);
phat2   = Xhat2(1:3,5);
vP2     = reshapeT(P2);
vxhat2  = [vRhat2;phat2;vhat2;vP2];


 


% output
tauplus = 1; % 1 second
xplus=[vx;vxhat1;vxhat2;tauplus;time];

end