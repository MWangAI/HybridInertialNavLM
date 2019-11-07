clear all
close all

global pointsNum g e3 Kn pointsI kR w bw BiasEnable
% System Parameters
g      = 9.8;
e3     = [0 0 1]';
bw     = [-0.1 0.02 0.02]';
% Set BiasEnable=1 to enable biase estimaiton and otherwise BiasEnable=0
BiasEnable = 1; 
% inertial frame landmarks   
pointsNum = 4;
Kn        = eye(pointsNum)/pointsNum;
% pointsI   = [-1.0258    1.4452    0.6591   -0.0693;
%              -2.4512   -0.8697    0.3664   -1.4358;
%               2.4979    0.2354   -2.6285    0.0542]; 
pointsI   = 2*randn(3,pointsNum);%[randn(2,pointsNum); zeros(1,pointsNum)];           
% pointsC   = sum(pointsI*Kn,2)/pointsNum; % center of landmarks 

kc    = sum(diag(Kn));
Temp  = pointsI-sum(pointsI*Kn,2)/sum(diag(Kn)); % p-pc
M     = Temp * Kn * Temp';

[U,E]  = eig(M);
Mbar   = (trace(M)*eye(3)-M)/2;
kR     = 5;

 


 
%%  initial conditions
w     = 0.6;
u     = e3;
R0    = expm(0*Skew(e3));
p0    = 10*[1 0 1]';
v0    = 10*w*[0 1 0]';
vR0   = reshape(R0,9,1);
x0    = [vR0;p0;v0];
 

% Initialization for observer
Rq     = expm(0.1*pi*Skew(U(:,3)));
Rhat0  = Rq'*R0;
bwhat0 = zeros(3,1);
vhat0  = 0*v0;
phat0  = 0*p0 + 0.0*randn(3,1);
vRhat0 = reshape(Rhat0,9,1);
xhat0  = [vRhat0;phat0;vhat0;bwhat0]; 





% HINO-CRE
P0    = 0.5*eye(6);
vP0   = reshapeT(P0);  


xIn = [x0;xhat0;xhat0;vP0];   
TSPAN=[0 15];

[Tout,Yout] = ode45(@myfun,TSPAN,xIn);

Rout1     = Yout(:,1:9);
pout1     = Yout(:,10:12);
vout1     = Yout(:,13:15);  
index     = 15; % HINO
Rhatout1  = Yout(:,index+1:index+9);
phatout1  = Yout(:,index+10:index+12);
vhatout1  = Yout(:,index+13:index+15);
bwhatout1 = Yout(:,index+16:index+18);
index     = 33;  % HINO-CRE
Rhatout2  = Yout(:,index+1:index+9);
phatout2  = Yout(:,index+10:index+12);
vhatout2  = Yout(:,index+13:index+15);
bwhatout2 = Yout(:,index+16:index+18);
% lP1      = length(A1)*(length(A1)+1)/2;



figure
% plot3(r(1,:),r(2,:),r(3,:),'k*'), hold on
plot3(pout1(:,1),pout1(:,2),pout1(:,3),'y--','linewidth',1), hold on 
% HINO
plot3(phatout1(:,1),phatout1(:,2),phatout1(:,3),'b-','linewidth',1)
% HINO-CRE
plot3(phatout2(:,1),phatout2(:,2),phatout2(:,3),'r-','linewidth',1)
grid on
% legend('true trajectory','A-HINO', 'iEKF')
% % axis([-11 11 -11 11 0 6])


perror = zeros(2,length(Tout));
Rerror = zeros(2,length(Tout));
berror = zeros(2,length(Tout));
for i=1:length(Tout)
    vR     = Rout1(i,:);
    R      = reshape(vR,3,3);
    p      = pout1(i,:)';
    v      = vout1(i,:)';
%     X      = [R v p;zeros(2,3) eye(2)];
    
    vRhat1 = Rhatout1(i,:);
    Rhat1  = reshape(vRhat1,3,3);
    phat1  = phatout1(i,:)';
    vhat1  = vhatout1(i,:)'; 
    bwhat1 = bwhatout1(i,:)';
%     XhatInv1  = [Rhat1' -Rhat1'*vhat1 -Rhat1'*phat1;zeros(2,3) eye(2)];
    
    
    vRhat2 = Rhatout2(i,:);     
    Rhat2  = reshape(vRhat2,3,3);
    phat2  = phatout2(i,:)';
    vhat2  = vhatout2(i,:)';
    bwhat2 = bwhatout2(i,:)';
%     XhatInv2  = [Rhat2' -Rhat2'*vhat2 -Rhat2'*phat2;zeros(2,3) eye(2)];
    
%     verror1(:,i) = norm(v - vhat);
    perror(1,i) = norm(p - phat1);
    Rerror(1,i) = trace((eye(3)-R*Rhat1'))/4;%trace((eye(5)-X*XhatInv1)*(eye(5)-X*XhatInv1)');%
    berror(1,i) = norm(bwhat1-bw);
    perror(2,i) = norm(p - phat2);
    Rerror(2,i) = trace((eye(3)-R*Rhat2'))/4;
    berror(2,i) = norm(bwhat2-bw);

end

figure
subplot(3,1,1)
plot(Tout, Rerror(1,:),'b-',Tout, Rerror(2,:),'r-','linewidth',2);
% xlabel('$t$','interpreter','latex')
ylabel('$\|\tilde{R}\|_I^2$','interpreter','latex')
legend('HINO','HINO-CRE')
grid on
subplot(3,1,2)
plot(Tout, perror(1,:),'b-',Tout, perror(2,:),'r-','linewidth',2);
xlabel('$t$','interpreter','latex')
ylabel('$\|p-\hat{p}\|$','interpreter','latex')
grid on
subplot(3,1,3)
plot(Tout, berror(1,:),'b-',Tout, berror(2,:),'r-','linewidth',2);
xlabel('$t$','interpreter','latex')
ylabel('$\|b_\omega-\hat{b}_\omega\|$','interpreter','latex')
grid on


function xdot = myfun(t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  pointsNum g e3 Kn pointsI kR bw BiasEnable

%% state
% xIn = [x0;xhat0;vP10;xhat0;vP20;tau;time];  % 15+(15+3)+(15+45) +2
% real system states
vR     = x(1:9,1);
p      = x(10:12,1);
v      = x(13:15,1);
index  = 15;
vRhat1 = x(index+1:index+9,1);
phat1  = x(index+10:index+12,1);
vhat1  = x(index+13:index+15,1);
bwhat1 = x(index+16:index+18,1);
index  = 33;
vRhat2 = x(index+1:index+9,1);
phat2  = x(index+10:index+12,1);
vhat2  = x(index+13:index+15,1);
bwhat2 = x(index+16:index+18,1);
vP    = x(index+19:end,1);

 
% omega and accel measurements
R         = reshape(vR,3,3);
[a,omega] = imu(R,t);
a_y       = a + 0*randn(3,1);
omega_y   = omega + BiasEnable*bw+0*randn(3,1);
pointsB   = R'*(pointsI - p) + 0*randn(size(pointsI));
r         = [pointsI;zeros(1,pointsNum);ones(1,pointsNum)];
b         = [pointsB;zeros(1,pointsNum);ones(1,pointsNum)];

pC = sum(pointsI*Kn,2)/pointsNum; % center of landmarks 
Xc      = [eye(3) zeros(3,1) pC;zeros(2,3) eye(2)];
XcInv   = [eye(3) zeros(3,1) -pC;zeros(2,3) eye(2)];

% dynamics of real system 
DR    = R*Skew(omega);
Dp    = v;
Dv    = -g*e3 + R*a;
DvR   = reshape(DR,9,1);
Dx    = [DvR;Dp;Dv];

% HINO
kp = 3;
kv = 3;
kw   = 1;
Rhat1   = reshape(vRhat1,3,3);
Xhat1   = [Rhat1 vhat1 phat1;zeros(2,3) eye(2)];
% K       = [kR*eye(3) zeros(3,2);zeros(1,5);zeros(1,3) kv kp];
% Temp    = XcInv*(r-Xhat1*b)*Kn*r'*XcInv'*K;
% Temp = [kR*(Temp(1:3,1:3)-Temp(1:3,1:3)')/2 kv*Temp(1:3,5) kp*Temp(1:3,5);zeros(2,5)];
Temp    = GainMap(XcInv*(r-Xhat1*b)*Kn*r'*XcInv',kR,kv,kp);
Delta   = Xc*Temp*XcInv;

DRhat1  = Rhat1*Skew(omega_y-BiasEnable*bwhat1)+Delta(1:3,1:3)*Rhat1;
Dphat1  = vhat1+Delta(1:3,1:3)*(phat1-pC) + Delta(1:3,5);
Dvhat1  = -g*e3 + Rhat1*a_y + Delta(1:3,1:3)*vhat1 + Delta(1:3,4);
psiR =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
Dbwhat1 = BiasEnable*(-kw*Rhat1'*psiR);

DvRhat1 = reshape(DRhat1,9,1);
Dxhat1  = [DvRhat1;Dphat1;Dvhat1;Dbwhat1];

% HINO-CRE
A = [-Skew(omega_y-BiasEnable*bwhat2) eye(3);zeros(3) -Skew(omega_y-BiasEnable*bwhat2)];
C = [eye(3) zeros(3)];

P      = reshapeT(vP);
Rhat2   = reshape(vRhat2,3,3);
Xhat2   = [Rhat2 vhat2 phat2;zeros(2,3) eye(2)];
V      = 1*diag([1*ones(1,3) 1*ones(1,3)]);
Q      = 1.0e-1*eye(3);
KP     = 1*P*C'/Q;
Kp     = Rhat2*KP(1:3,1:3)*Rhat2'/sum(diag(Kn));%Rhat2*KP(1:3,1:3)*Rhat2';%
Kv     = Rhat2*KP(4:6,1:3)*Rhat2'/sum(diag(Kn)); %Rhat2*KP(4:6,1:3)*Rhat2';%
% Temp   = XcInv*(r-Xhat2*b)*Kn*r'*XcInv';
% Temp    = [kR*Temp(1:3,1:3) Kv*Temp(1:3,5) Kp*Temp(1:3,5);zeros(2,5)];
Temp   = GainMap(XcInv*(r-Xhat2*b)*Kn*r'*XcInv',kR,Kv,Kp);
Delta   = Xc*Temp*XcInv;

DRhat2  = Rhat2*Skew(omega_y-BiasEnable*bwhat2)+Delta(1:3,1:3)*Rhat2;
Dphat2  = vhat2+Delta(1:3,1:3)*(phat2-pC) + Delta(1:3,5);
Dvhat2  = -g*e3 + Rhat2*a_y + Delta(1:3,1:3)*vhat2 + Delta(1:3,4);
psiR =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
Dbwhat2 = BiasEnable*(-kw*Rhat2'*psiR);
DP     = A*P+P*A' + V - (P*C'/Q)*C*P;

DvRhat2 = reshape(DRhat2,9,1);
DvP    = reshapeT(DP);
Dxhat2  = [DvRhat2;Dphat2;Dvhat2;Dbwhat2;DvP];

% output
xdot=[Dx;Dxhat1;Dxhat2];

end