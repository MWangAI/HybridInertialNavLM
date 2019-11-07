function [Tout,error,phatout]=HINO(TSPAN,xIn,options)
% HINO
global bw
JSPAN = [0 30];
% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;
% options = odeset('RelTol',1e-1,'MaxStep',.01);
[Tout,J,Yout] = HyEQsolver(@f1,@g1,@C1,@D1,xIn,TSPAN,JSPAN,rule,options); 


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


function [value] = C1(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                
%
% Description: Flow set
% Return 0 if outside of C, and 1 if inside C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   Kn pointsI
%% state 

Q     = x(1:4,1);
p      = x(5:7,1);
% v      = x(8:10,1);
index  = 10;
Qhat  = x(index+1:index+4,1);
% phat  = x(index+5:index+7,1);
% vhat  = x(index+8:index+10,1);
% bwhat = x(index+11:index+13,1);
% t     = x(end,1); 
R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat');

kc        = sum(diag(Kn));
pointsB   = R'*(pointsI - p) + 0*randn(size(pointsI));
pointsNum = size(pointsI,2);
pC        = sum(pointsI*Kn,2)/kc;
yC        = sum(pointsB*Kn,2)/kc;
error     = (pointsI - pC)-Rhat*(pointsB-yC);
Phi       = 0.5*trace(error'*error*Kn);

Temp  = pointsI-pC; % p-pc
M     = zeros(3);
for i=1:pointsNum
   M = M + Kn(i,i)*Temp(:,i)*Temp(:,i)';
end
[U,E]  = eig(M);
lambda1 = min(diag(E));
lambda3 = max(diag(E));
lambda2 = trace(E)-lambda1-lambda3;
if lambda1==lambda3 
    DeltaM = lambda1*2/3;
elseif (lambda1==lambda2)||(lambda2==lambda3)
    DeltaM = min(2*lambda2,trace(E)-2*lambda2);
else
    DeltaM = lambda1 + lambda2;
end

theta = 0.8*pi;
delta = 0.3*(1-cos(theta))*DeltaM;
[~,m] = size(U);
Phiq  = zeros(1,m);
for i=1:m
   Rq      = expm(theta*Skew(U(:,i)));
%    pq      = (eye(3)-Rq)*pC;
%    vq      = zeros(3,1);
%    XqInv   = [Rq' -Rq'*vq -Rq'*pq;zeros(2,3) eye(2)];
   error   = (pointsI - pC)-Rq'*Rhat*(pointsB-yC);
   Phiq(i) = 0.5*trace(error'*error*Kn); 
end 


if (Phi-min(Phiq)<=delta)
    value = 1;
else 
    value = 0;
end

end


function [inside] = D1(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                
%
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   Kn pointsI LandMnoise
%% state 
Q     = x(1:4,1);
p      = x(5:7,1);
% v      = x(8:10,1);
index  = 10;
Qhat  = x(index+1:index+4,1);
% phat  = x(index+5:index+7,1);
% vhat  = x(index+8:index+10,1);
% bwhat = x(index+11:index+13,1);
% t     = x(end,1); 
R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat');


kc        = sum(diag(Kn));
pointsB   = R'*(pointsI - p) + LandMnoise*randn(size(pointsI));
pointsNum = size(pointsI,2);
pC        = sum(pointsI*Kn,2)/kc;
yC        = sum(pointsB*Kn,2)/kc;
error     = (pointsI - pC)-Rhat*(pointsB-yC);
Phi       = 0.5*trace(error'*error*Kn);

Temp  = pointsI-pC; % p-pc
M     = zeros(3);
for i=1:pointsNum
   M = M + Kn(i,i)*Temp(:,i)*Temp(:,i)';
end
[U,E]  = eig(M);
lambda1 = min(diag(E));
lambda3 = max(diag(E));
lambda2 = trace(E)-lambda1-lambda3;
if lambda1==lambda3 
    DeltaM = lambda1*2/3;
elseif (lambda1==lambda2)||(lambda2==lambda3)
    DeltaM = min(2*lambda2,trace(E)-2*lambda2);
else
    DeltaM = lambda1 + lambda2;
end

theta = 0.8*pi;
delta = 0.3*(1-cos(theta))*DeltaM;
[~,m] = size(U);
Phiq  = zeros(1,m);
for i=1:m
   Rq      = expm(theta*Skew(U(:,i)));
%    pq      = (eye(3)-Rq)*pC;
%    vq      = zeros(3,1);
%    XqInv   = [Rq' -Rq'*vq -Rq'*pq;zeros(2,3) eye(2)];
   error   = (pointsI - pC)-Rq'*Rhat*(pointsB-yC);
   Phiq(i) = 0.5*trace(error'*error*Kn); 
end 


if (Phi-min(Phiq)>=delta)
    inside = 1;
else 
    inside = 0;
end

end


function xdot=f1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   g e3 Kn pointsI kR bw BiasEnable kw kv kp LandMnoise IMUnoise
% real system states
Q     = x(1:4,1);
p      = x(5:7,1);
v      = x(8:10,1);
index  = 10;
Qhat  = x(index+1:index+4,1);
phat  = x(index+5:index+7,1);
vhat  = x(index+8:index+10,1);
bwhat = x(index+11:index+13,1);
t     = x(end,1); 
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
% DR    = R*Skew(omega);
Dp    = v;
Dv    = -g*e3 + R*a;
Dx    = [dQ';Dp;Dv];

% HINO
% kv   = 25;
% kp   = 25;
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
xdot=[Dx;Dxhat;1];
end

function  xplus = g1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   Kn pointsI LandMnoise
%% state 
Q     = x(1:4,1);
p      = x(5:7,1);
v      = x(8:10,1);
index  = 10;
Qhat  = x(index+1:index+4,1);
phat  = x(index+5:index+7,1);
vhat  = x(index+8:index+10,1);
bwhat = x(index+11:index+13,1);
t     = x(end,1); 
R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat');

kc        = sum(diag(Kn)); 
pointsB   = R'*(pointsI - p) + LandMnoise*randn(size(pointsI));
pointsNum = size(pointsI,2);
pC        = sum(pointsI*Kn,2)/kc;
yC        = sum(pointsB*Kn,2)/kc;
% error     = (pointsI - pC)-Rhat*(pointsB-yC);
% Phi       = 0.5*trace(error'*error*Kn);

Temp  = pointsI-pC; % p-pc
M     = zeros(3);
for i=1:pointsNum
   M = M + Kn(i,i)*Temp(:,i)*Temp(:,i)';
end
[U,E]  = eig(M);
% lambda1 = min(diag(E));
% lambda3 = max(diag(E));
% lambda2 = trace(E)-lambda1-lambda3;
% if lambda1==lambda3 
%     DeltaM = lambda1*2/3;
% elseif (lambda1==lambda2)||(lambda2==lambda3)
%     DeltaM = min(2*lambda2,trace(E)-2*lambda2);
% else
%     DeltaM = lambda1 + lambda2;
% end

theta = 0.8*pi;
% delta = 0.3*(1-cos(theta))*DeltaM;
[~,m] = size(U);
Phiq  = zeros(1,m);
for i=1:m
   Rq      = expm(theta*Skew(U(:,i)));
%    pq      = (eye(3)-Rq)*pC;
%    vq      = zeros(3,1);
%    XqInv   = [Rq' -Rq'*vq -Rq'*pq;zeros(2,3) eye(2)];
   error   = (pointsI - pC)-Rq'*Rhat*(pointsB-yC);
   Phiq(i) = 0.5*trace(error'*error*Kn); 
end 

[~,index]  = min(Phiq);

Rq        = expm(theta*Skew(U(:,index)));
pq        = (eye(3)-Rq)*pC;
vq        = zeros(3,1);
XqInv     = [Rq' -Rq'*vq -Rq'*pq;zeros(2,3) eye(2)];
Xhat      = [Rhat vhat phat;zeros(2,3) eye(2)];
Xhat_plus = XqInv*Xhat;

Qq = [cos(theta/2) sin(theta/2)*U(:,index)'];
QqInv = quatinv(Qq);
Qhat_plus = quatmultiply(QqInv,Qhat');
% Rhat_plus  = Xhat_plus(1:3,1:3);
% vRhat_plus = reshape(Rhat_plus,9,1);
vhat_plus  = Xhat_plus(1:3,4);
phat_plus  = Xhat_plus(1:3,5);

xplus = [Q;p;v;Qhat_plus';phat_plus;vhat_plus;bwhat;t];

end