clear
close all
clc
 
global  pointsNum g e3 Kn pointsI bw ba kR kw kv kp w IMUnoise LandMnoise
% System Parameters
g      = 9.8;
e3     = [0 0 1]';
bw     = [-0.1 0.02 0.02]';%[-0.002229	0.0207	0.07635]';
ba     = [-0.012492	0.547666 0.069073]';

% inertial frame landmarks    
pointsNum = 6;
Kn        = eye(pointsNum)/pointsNum;
% pointsI   = [-1.0258    1.4452    0.6591   -0.0693;
%              -2.4512   -0.8697    0.3664   -1.4358;
%               2.4979    0.2354   -2.6285    0.0542]; 
pointsI  = [0.9564   -5.0507    1.6929   -3.9198   -1.5798    2.7664;
            0.9907    3.6939   -0.0476    1.4497   -0.8332   -0.8695;
            4.3994    1.1337    1.4279    1.3990    1.1610    0.6336];
% pointsI   = 2*randn(3,pointsNum);%[randn(2,pointsNum); zeros(1,pointsNum)];           
% pointsC   = sum(pointsI*Kn,2)/sum(diag(Kn)); % center of landmarks 

kc    = sum(diag(Kn));
Temp  = pointsI-sum(pointsI*Kn,2)/pointsNum; % p-pc
M     = Temp * Kn * Temp';

[U,E]  = eig(M);
Mbar   = (trace(M)*eye(3)-M)/2;
kR = 2;
 
 
IMUnoise = [0.1 0.1];
LandMnoise = 0.1;

 
 

 
%%  initial conditions
w     = 0.5;
u     = e3;
Q0    = [1 0 0 0];
p0    = 10*[1 0 1]';
v0    = 10*w*[0 1 0]';
x0    = [Q0';p0;v0];
t     = 0;

% Initialization for observer
theta  = 0.99*pi;
vhat0  = 0*v0+ 0.0*randn(3,1);
phat0  = 0*p0 + 0.0*randn(3,1);
 
Qhat0  = [cos(theta/2) sin(theta/2)*U(:,3)'];
bwhat0 = zeros(3,1); 
xhat0  = [Qhat0';phat0;vhat0;bwhat0]; 
bahat0 = zeros(3,1);


% HINO-CRE
P0    = 1*diag([1*ones(1,3) 1*ones(1,3) 1*ones(1,3)]);
vP0   = reshapeT(P0);  

%% 
 
TSPAN=[0 30];
options = odeset('RelTol',1e-1,'MaxStep',.01);



xIn3 = [x0;xhat0;bahat0;vP0;t];  
[Tout1,error1,phatout1]=HINOCRE(TSPAN,xIn3,options);
 
 
 
    
figure
pout = 10*[cos(w*Tout1) sin(w*Tout1) ones(size(Tout1))];
plot3(pout(:,1),pout(:,2),pout(:,3),'y--','linewidth',1), hold on
plot3(phatout1(:,1),phatout1(:,2),phatout1(:,3),'b-','linewidth',1)
grid on
%     axis([-11 11 -11 11 0 16])

legend('True trajectory','Hybrid-CRE2')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','3dTrajectory2.eps')


figure
subplot(2,2,1)
plot(Tout1, error1(1,:),'b-','linewidth',2);
xlabel('$t(s)$','interpreter','latex') 
ylabel('$|\tilde{R}|_I$','interpreter','latex')
    legend('Hybrid-CRE2')
grid on 
subplot(2,2,2)
plot(Tout1, error1(2,:),'b-','linewidth',2);
% xlabel('$t$','interpreter','latex')
ylabel('$\|p-\hat{p}\|$','interpreter','latex')
grid on
xlabel('$t(s)$','interpreter','latex')
subplot(2,2,3)
plot(Tout1, error1(3,:),'b-','linewidth',2);
xlabel('$t(s)$','interpreter','latex')
ylabel('$\|v-\hat{v}\|$','interpreter','latex')
grid on
xlabel('$t(s)$','interpreter','latex')
subplot(2,2,4)
yyaxis left
plot(Tout1, error1(4,:),'b-','linewidth',2), grid on
ylabel('$\|b_\omega-\hat{b}_\omega\|$','interpreter','latex')
yyaxis right
plot(Tout1, error1(5,:),'r-','linewidth',2), grid on
ylabel('$\|b_a-\hat{b}_a\|$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex') 
 


 


 