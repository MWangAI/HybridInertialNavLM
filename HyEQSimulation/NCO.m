function [Tout,error,pout,phatout]=NCO(TSPAN,xIn,options)
% Nonlinear continuous observer
global bw

% options = odeset('RelTol',1e-1,'MaxStep',.01);
[Tout,Yout] = ode45(@myfun,TSPAN,xIn,options);

Qout     = Yout(:,1:4);
pout     = Yout(:,5:7);
vout     = Yout(:,8:10);  
index     = 10; % HINO
Qhatout  = Yout(:,index+1:index+4);
phatout  = Yout(:,index+5:index+7);
vhatout  = Yout(:,index+8:index+10);
bwhatout = Yout(:,index+11:index+13);



error = zeros(4,length(Tout)); 
for i=1:length(Tout)
    Q     = Qout(i,:);
    R      = quat2rotm(Q);
    p      = pout(i,:)';
    v      = vout(i,:)';
%     X      = [R v p;zeros(2,3) eye(2)];
    
    Qhat = Qhatout(i,:);
    Rhat  = quat2rotm(Qhat);
    phat  = phatout(i,:)';
    vhat  = vhatout(i,:)'; 
    bwhat = bwhatout(i,:)';
%     XhatInv1  = [Rhat1' -Rhat1'*vhat1 -Rhat1'*phat1;zeros(2,3) eye(2)];
    
  
    error(3,i) = norm(v - vhat);
    error(2,i) = norm(p - phat);
    error(1,i) = (trace((eye(3)-R*Rhat'))/4);%trace((eye(5)-X*XhatInv1)*(eye(5)-X*XhatInv1)');%
    error(4,i) = norm(bwhat-bw);    

end

end

% clearvars  Yout Rout vout pout Rhatout phatout vhatout bwhatout


function xdot = myfun(t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global   g e3 Kn pointsI kR bw BiasEnable kw kv kp LandMnoise IMUnoise

%% state 
% real system states

Q     = x(1:4,1);
p      = x(5:7,1);
v      = x(8:10,1);
index  = 10;
Qhat  = x(index+1:index+4,1);
phat  = x(index+5:index+7,1);
vhat  = x(index+8:index+10,1);
bwhat = x(index+11:index+13,1);
R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat');
 

 
% omega and accel measurements
 
[a,omega] = imu(R,t);
a_y       = a + IMUnoise(1)*randn(3,1);
omega_y   = omega + BiasEnable*bw+IMUnoise(2)*randn(3,1);
pointsB   = R'*(pointsI - p) + LandMnoise*randn(size(pointsI));
pointsNum = size(pointsI,2);
r         = [pointsI;zeros(1,pointsNum);ones(1,pointsNum)];
b         = [pointsB;zeros(1,pointsNum);ones(1,pointsNum)];

pC = sum(pointsI*Kn,2)/sum(diag(Kn)); % center of landmarks 
Xc      = [eye(3) zeros(3,1) pC;zeros(2,3) eye(2)];
XcInv   = [eye(3) zeros(3,1) -pC;zeros(2,3) eye(2)];

% dynamics of real system 
Qwhat = [0;omega];
dQ = 0.5*quatmultiply(Q',Qwhat');
Dp    = v;
Dv    = -g*e3 + R*a;
Dx    = [dQ';Dp;Dv];

% NCO
% kv   = 20;
% kp   = 20;
 
Xhat   = [Rhat vhat phat;zeros(2,3) eye(2)];
% K       = [kR*eye(3) zeros(3,2);zeros(1,5);zeros(1,3) kv kp];
% Temp    = XcInv*(r-Xhat1*b)*Kn*r'*XcInv'*K;
% Temp = [kR*(Temp(1:3,1:3)-Temp(1:3,1:3)')/2 kv*Temp(1:3,5) kp*Temp(1:3,5);zeros(2,5)];
Temp    = GainMap(XcInv*(r-Xhat*b)*Kn*r'*XcInv',kR,kv,kp);
Delta   = Xc*Temp*XcInv;
psiR =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
% DRhat  = Rhat*Skew(omega_y-BiasEnable*bwhat)+Delta(1:3,1:3)*Rhat;
Qwhat = [0;omega_y-BiasEnable*bwhat+ Rhat'*psiR];
dQhat = 0.5*quatmultiply(Qhat',Qwhat');
Dphat  = vhat+Delta(1:3,1:3)*(phat-pC) + Delta(1:3,5);
Dvhat  = -g*e3 + Rhat*a_y + Delta(1:3,1:3)*vhat + Delta(1:3,4);

Dbwhat = BiasEnable*(-kw*Rhat'*psiR);
Dxhat  = [dQhat';Dphat;Dvhat;Dbwhat];

 

% output
xdot=[Dx;Dxhat];

end