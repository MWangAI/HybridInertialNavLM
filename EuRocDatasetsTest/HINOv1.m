clear all
close all
clc

% Navigation on SE_2(3) with biased angular velocity
% Limited number of featrues in the pointTracker!
% for example 30. This method saves time for the feature tracking
%% Read dataset

% Load stereo camera intrinsic parameters (calibrated using MATLAB)
load('stereoParamsdataset.mat'); 

% Set the path of different experiments
setpath = 1;
switch  setpath
    case 1
        path        = 'V1_01_easy'; 
    case 2
        path        = 'V1_02_medium'; 
    case 3
        path        = 'V2_01_easy'; 
    case 4
        path        = 'V2_02_medium'; 
end

% Read the dataset
[IMURead, ViconRead, Groundtruth, CamRead]= dataRead(path); 
bw        = mean(Groundtruth(:,12:14)); %averaged gyro bias from groundtruth
ba        = mean(Groundtruth(end-100:end,15:17)); %averaged accel bias from groundtruth
TimerIMU  = IMURead(:,1)*1e-9;  %nanoseconds to seconds 
TimerCam  = CamRead(:,1)*1e-9;  
% transfromation from body frame to cam0 frame 
TB2C0 = [0.0148655429818, -0.999880929698, 0.00414029679422, -0.0216401454975;
          0.999557249008, 0.0149672133247, 0.025715529948, -0.064676986768;
          -0.0257744366974, 0.00375618835797, 0.999660727178, 0.00981073058949;
          0.0, 0.0, 0.0, 1.0];
 
%% Keyframe initialization
numFeatures = 60;
index       = 1;
[imageLeft,imageRight]                       = stereoRead(path,CamRead,index,stereoParams);
% pointTracker = vision.PointTracker;
pointTracker = vision.PointTracker('BlockSize',[11 11],'MaxBidirectionalError',1); %Reduce the blocksize can improve the accuracy
[pointsKeyR,pointsKeyL,pointTracker] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker);
% Show matched features
% figure; ax = axes;
% showMatchedFeatures(imageLeft,imageRight,pointsKeyL,pointsKeyR,'blend','Parent',ax);
% % title(ax, 'Candidate point matches');
% legend(ax, 'Matched points left','Matched points right');
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','Features.eps')
%% Initialize observer

g  = 9.81;
e3 = [0,0,1]';
 
P  = 0.1*eye(6); 
kR = 0.15;
kv = 0.5;
kp = 0.5;
kw = 30;
Q  = 20*eye(3);
V  = 10*diag([1*ones(1,3) 1*ones(1,3)]); 
    
    
p  = Groundtruth(1,2:4)';
R  = quat2rotm(Groundtruth(1,5:8));
v  = Groundtruth(1,9:11)';
 



Rkey = R;
pkey = p;
vkey = v; 
points3D  = funPoints3D(pointsKeyL,pointsKeyR,stereoParams,TB2C0);
pointsKey = Rkey*points3D+pkey;  %inertial frame landmarks  

u      = e3;%randn(3,1);
theta  = 0.99*pi;%(0.5-rand(1))*0.4*pi;
Rq     = expm(Skew(u)*theta);
Xhat   = [Rq*R 0*v+0.0*randn(3,1) 0*p+0.0*randn(3,1);
         zeros(2,3) eye(2)];
% Xhat   = [Rq*R 1*v+0.0*randn(3,1) 1*p+0.0*randn(3,1);
%          zeros(2,3) eye(2)];
% Xhat   = [eye(3) 0*v 0*p;zeros(2,3) eye(2)];

bwhat  = zeros(3,1);
Xhat2  = Xhat;
bwhat2 = bwhat;


%% Loop
Tfinal     = ceil(length(TimerIMU)*0.96);%20000;%
Pestimate  = zeros(3,Tfinal);
Pestimate2 = zeros(3,Tfinal); 
normR      = zeros(2,Tfinal);
normV      = zeros(2,Tfinal); 
normbw     = zeros(2,Tfinal);
Numlandmark    = zeros(size(CamRead));
keytime        = zeros(size(CamRead));
keytime(index) = 1;
for k=1:1:Tfinal
    if mod(floor(100*k/Tfinal),5)==0
    clc
    disp(['Complete  ' num2str(floor(100*k/Tfinal)) '%'])
    end
    % errors
    p       = Groundtruth(k,2:4)';
    R       = quat2rotm(Groundtruth(k,5:8));
    v       = Groundtruth(k,9:11)';     
    normR(1,k)      = trace(eye(3)-R*Xhat(1:3,1:3)')/4; 
    normR(2,k)      = trace(eye(3)-R*Xhat2(1:3,1:3)')/4;
    normV(1,k)      = norm(v-Xhat(1:3,4)); 
    normV(2,k)      = norm(v-Xhat2(1:3,4));
    normbw(1,k)     = norm(bw'-bwhat);
    normbw(2,k)     = norm(bw'-bwhat2);
    Pestimate(:,k)  = Xhat(1:3,5)'; 
    Pestimate2(:,k) = Xhat2(1:3,5)';
%         propagation        
    dT    = TimerIMU(k+1)-TimerIMU(k);            
    omega = IMURead(k,2:4)' ;
    a     = IMURead(k,5:7)' ;        
    
%     Adaptive gains
    Rhat = Xhat(1:3,1:3);
    vhat = Xhat(1:3,4); 
    phat = Xhat(1:3,5);     
    Rhat = Rhat*expm(dT*Skew(omega-bwhat));
    phat = phat + vhat*dT;
    vhat = vhat + (-g*e3+Rhat*(a-1*ba'))*dT;
    Xhat = [Rhat vhat phat;zeros(2,3) eye(2)];
 
    
%     Constant gains
    Rhat2 = Xhat2(1:3,1:3);
    vhat2 = Xhat2(1:3,4);
    phat2 = Xhat2(1:3,5);
    Rhat2 = Rhat2*expm(dT*Skew(omega-bwhat2));
    phat2 = phat2 + vhat2*dT;
    vhat2 = vhat2 + (-g*e3+Rhat2*(a-1*ba'))*dT;
    Xhat2 = [Rhat2 vhat2 phat2;zeros(2,3) eye(2)]; 
 
%    
    
    A = [-Skew(omega-bwhat) eye(3);zeros(3) -Skew(omega-bwhat)];
    C = [eye(3) zeros(3)];
    
    Ad    = expm(A*dT);
    P     = Ad*P*Ad'+V; %P + (A*P*A' + V)*dT;%
    
    
    
    
    if TimerIMU(k+1)>= TimerCam(index+1,1)
        index = index +1;               
        [imageLeft,imageRight]             = stereoRead(path,CamRead,index,stereoParams);
        [pointsNewL,pointsNewR,isFoundand] = stereoPointTrack(imageLeft,imageRight,pointTracker);        
        % tracked features locations
        pointsNewR = pointsNewR(isFoundand,:);
        pointsNewL = pointsNewL(isFoundand,:);
%         pointsNewL = pointsNewL + 5*randn(size(pointsNewL)); %add noise for analysis
%         pointsNewR = pointsNewR + 5*randn(size(pointsNewR));
        
        pointsI    = pointsKey(:,isFoundand);
        pointsB   = funPoints3D(pointsNewL,pointsNewR,stereoParams,TB2C0);  
        [pointsI,pointsB,pointsNum] = pointsSelection(pointsI,pointsB); 
        Numlandmark(index)= size(pointsI,2);
%         pointsNum = size(pointsI,2);
        
        if (pointsNum>=3)
            % 
            Kn      = eye(pointsNum)/pointsNum;
            pointsC = sum(pointsI*Kn,2); % center of landmarks
            M       = (pointsI-pointsC)*Kn*(pointsI-pointsC)';
            Mbar   = (trace(M)*eye(3)-M)/2;
            if max(eig(Mbar))>0.01
                kR     = 0.4/max(eig(Mbar));
            end
            
            
            r = [pointsI;zeros(1,pointsNum);ones(1,pointsNum)];
            b = [pointsB;zeros(1,pointsNum);ones(1,pointsNum)];
                
            Xc      = [eye(3) zeros(3,1) pointsC;zeros(2,3) eye(2)];
            XcInv   = [eye(3) zeros(3,1) -pointsC;zeros(2,3) eye(2)]; 
             
            
            % Adaptive gain    
			[XqInv,flag] = RestX(Xhat,r,b,Kn);
            if  flag == 1
                    Xhat = XqInv*Xhat;
                    
            end
            
            KP     = P*C'/(C*P*C'+Q);
            Kp     = Rhat*KP(1:3,1:3)*Rhat'/sum(diag(Kn));%Rhat2*KP(1:3,1:3)*Rhat2';%
            Kv     = Rhat*KP(4:6,1:3)*Rhat'/sum(diag(Kn)); %Rhat2*KP(4:6,1:3)*Rhat2';%             
            Temp    = GainMap(XcInv*(r-Xhat*b)*Kn*r'*XcInv',kR,Kv,Kp); 
            Delta   = Xc*Temp*XcInv;
            
%             Rhat    = expm(Delta(1:3,1:3))*Rhat;
%             phat    = phat + Delta(1:3,1:3)*(phat-pointsC) + Delta(1:3,5);
%             vhat    = vhat + Delta(1:3,1:3)*vhat + Delta(1:3,4);
%             Xhat    = [Rhat vhat phat;zeros(2,3) eye(2)];
            XhatInv = InvX(Xhat); 
            Xhat    = Xhat*expm(XhatInv*Delta*Xhat); % matrix form
            psiR   =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
            bwhat = bwhat + dT*(-kw*Rhat'*psiR);
            
            P       = (eye(size(P))-KP*C)*P;
 
            
            % Constant gain 
			[XqInv,flag] = RestX(Xhat2,r,b,Kn);
            if  flag == 1
                    Xhat = XqInv*Xhat2;
                    
            end
			
              
            Temp     = GainMap(XcInv*(r-Xhat2*b)*Kn*r'*XcInv',kR,kv,kp);                    
            Delta    = Xc*Temp*XcInv;            
%             Rhat2    = expm(Delta(1:3,1:3))*Rhat2;
%             phat2    = phat2 + Delta(1:3,1:3)*(phat2-pointsC) + Delta(1:3,5);
%             vhat2    = vhat2 + Delta(1:3,1:3)*vhat2 + Delta(1:3,4);
%             Xhat2    = [Rhat2 vhat2 phat2;zeros(2,3) eye(2)];
            XhatInv2 = InvX(Xhat2); 
            Xhat2    = Xhat2*expm(XhatInv2*Delta*Xhat2);
            psiR =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
            bwhat2 = bwhat2 + dT*(-kw*Rhat2'*psiR);  
            
        end %endif
        % Creat new Keyframe if pointsNum is less than 6     
        minFeatures = 6;
        if pointsNum <= minFeatures %size(find(isFoundand==1), 1) < minFeatures       
            [pointsKeyR,pointsKeyL,pointTracker] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker);                 
            pkey      = Groundtruth(k,2:4)';
            Rkey      = quat2rotm(Groundtruth(k,5:8));
            points3D  = funPoints3D(pointsKeyL,pointsKeyR,stereoParams,TB2C0);
            pointsKey = Rkey*points3D+pkey;  %inertial frame landmarks  
            keytime(index) = 1;

         end %endif
    end %endif
    
    
    
end %endfor
 
release(pointTracker);
 
clearvars -except TimerIMU Groundtruth Tfinal Pestimate Pestimate2 Numlandmark TimerCam keytime normR normbw normV
%% Plot
PGrdtruth = Groundtruth(1:Tfinal,2:4)';
errorP     = Pestimate - PGrdtruth;
errorPnorm = sqrt(sum(errorP.^2,1));
RMSE       = sqrt(mean(errorPnorm.^2)) 

errorP2     = Pestimate2 - PGrdtruth;
errorPnorm2 = sqrt(sum(errorP2.^2,1));
RMSE2       = sqrt(mean(errorPnorm2.^2)) 

Time = (TimerIMU(1:Tfinal)-TimerIMU(1));
figure 
subplot(2,2,1)
plot(Time,sqrt(normR(2,:))','r-',Time,sqrt(normR(1,:))','b-','linewidth',1), grid on
legend('HINO', 'HINO-CRE')
ylabel('$|\tilde{R}|_I$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
xlim([0 Time(end)])
% ylim([0 1.0])
%     axes('position',[.22 .72 .20 .22])
%     box on % put box around new pair of axes
%     indexOfInterest1  = (Time>=95) & (Time<=100); 
%     plot(Time(indexOfInterest1),sqrt(normR(2,indexOfInterest1))','r-',...
%         Time(indexOfInterest1),sqrt(normR(1,indexOfInterest1))','b-','linewidth',1), grid on
%     axis tight



subplot(2,2,2) 
plot(Time,errorPnorm2,'r-',Time,errorPnorm,'b-','linewidth',1), grid on
xlim([0 Time(end)])
% ylim([0 0.9]) 
ylabel('$\|p-\hat{p}\|$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
%     axes('position',[.72 .72 .20 .22])
%     box on % put box around new pair of axes
%     indexOfInterest1  = (Time>=95) & (Time<=100); 
%     plot(Time(indexOfInterest1), (errorPnorm2(indexOfInterest1))','r-',...
%         Time(indexOfInterest1), (errorPnorm(indexOfInterest1))','b-','linewidth',1), grid on
%     axis tight




subplot(2,2,3)
plot(Time,normV(2,:)','r-',Time,normV(1,:),'b-','linewidth',1), grid on
ylabel('$\|v-\hat{v}\|$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
xlim([0 Time(end)])
% ylim([0 1.6]) 




subplot(2,2,4)
plot(Time,normbw(2,:)','r-',Time,normbw(1,:),'b-','linewidth',1), grid on
ylabel('$\|b_\omega-\hat{b}_\omega\|$','interpreter','latex')
xlim([0 Time(end)])
xlabel('$t(s)$','interpreter','latex')
% ylim([0 0.16]) 

 set(gcf, 'Renderer', 'Painters');
 print('-depsc','SimulationError3.eps')
 
figure
plot(TimerCam-TimerIMU(1),Numlandmark,'.'), hold on
xlim([0 Time(end)])
Keyframetime = TimerCam(keytime==1);  
plot(Keyframetime-TimerIMU(1),zeros(size(Keyframetime)),'*')

figure 
plot3(PGrdtruth(1,:), PGrdtruth(2,:),PGrdtruth(3,:),'g--','linewidth',1), hold on
plot3(Pestimate2(1,:), Pestimate2(2,:),Pestimate2(3,:),'r-','linewidth',1) 
plot3(Pestimate(1,:), Pestimate(2,:),Pestimate(3,:),'b-','linewidth',1) 
% title(['RMSE=' num2str(RMSE)])
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
grid on
legend('True trajectory','HINO', 'HINO-CRE')
 set(gcf, 'Renderer', 'Painters');
 print('-depsc','3dTrajectory3.eps')

% set(gcf, 'Renderer', 'Painters');
% print('-depsc','E:\Dropbox (Personal)\Research Note\2-SE(3)\IEEE TAC 2018\3DTrajEx.eps')


% temp1 = sqrt(sum(Pestimate.^2,2));
% temp2 = sqrt(sum(Groundtruth(1:length(Pestimate),2:4).^2,2));
% RMSE = sqrt(mean((temp1-temp2).^2)); 

 

function [IMURead,ViconRead,Groundtruth,CamRead] = dataRead(path)

% IMU
fileID      = fopen([path '\imu0\data.csv']);
csvRawIMU   = textscan(fileID, '%f,%f,%f,%f,%f,%f,%f', 'headerLines', 1);
fclose(fileID);
IMURead     = cell2mat(csvRawIMU); 
% IMURead   = csvread([path '\imu0\data.csv'],1,0);


% Vicon
fileID      = fopen([path '\vicon0\data.csv']);
csvRawVicon = textscan(fileID, '%f,%f,%f,%f,%f,%f,%f,%f', 'headerLines', 1);
fclose(fileID);
ViconRead   = cell2mat(csvRawVicon);  
% ViconRead = csvread([path '\vicon0\data.csv'],1,0);

% state_groundtruth
fileID      = fopen([path '\state_groundtruth_estimate0\data.csv']);
csvRawState = textscan(fileID, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',...
                        'headerLines', 1);
fclose(fileID);
Groundtruth = cell2mat(csvRawState);  
% Groundtruth = csvread([path '\state_groundtruth_estimate0\data.csv'],1,0); 

fileID      = fopen([path '\cam0\data.csv']);
csvRawCam   = textscan(fileID, '%f,%s', 'headerLines', 1);
fclose(fileID);
CamRead     = [csvRawCam{1}];
% CamRead   = csvread([path '\cam0\data.csv'], '%d%d',1,1);


intTime     = intersect(IMURead(:,1),Groundtruth(:,1),'sorted');
timestart   = find(intTime(1,1)==IMURead(:,1));
timeend     = find(intTime(end,1)==IMURead(:,1));
IMURead     = IMURead(timestart:timeend,:);
timestart   = find(intTime(1,1)==Groundtruth(:,1));
timeend     = find(intTime(end,1)==Groundtruth(:,1));
Groundtruth = Groundtruth(timestart:timeend,:);

timestart   = find(CamRead(:,1)<Groundtruth(1,1), 1, 'last' );  %last max
timeend     = find(CamRead(:,1)>Groundtruth(end,1),1);          %first min
CamRead     = CamRead(timestart:timeend,:);

end

function [imageLeft,imageRight] = stereoRead(path,CamRead,k,stereoParams)
% read and rectify with stereoParams
imageLeft  = imread([path '\cam0\data\',num2str(CamRead(k),'%d'),'.png']);
imageRight = imread([path '\cam1\data\',num2str(CamRead(k),'%d'),'.png']);
[imageLeft,imageRight] = rectifyStereoImages(imageLeft,imageRight,stereoParams);
end

function [positionsB] = funPoints3D(pointsFeaturesL,pointsFeaturesR,stereoParams,TB2C0)
%     position on the frame of left cam 
positionCam0 = 0.001*triangulate(pointsFeaturesL,pointsFeaturesR,stereoParams);
positionCam0 = positionCam0';
%     positionCam0 = R'(positionsB-p)
R = TB2C0(1:3,1:3);
p = TB2C0(1:3,4);
positionsB = R*positionCam0 + p;    %R'*(positionCam0-p);%  
    
end

function [pointsKeyR,pointsKeyL,pointTracker] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker)

% feature detection on the right image of keyframe and set the points as
% pointTracker

points          = detectMinEigenFeatures(imageRight);
pointsKeyR      = points.Location;

release(pointTracker);
% Initialize the tracker with the initial point locations and the initial video frame.
initialize(pointTracker, pointsKeyR, imageRight);
% track features on the left image of keyframe
% [pointsKeyL, isFound] = step(pointTracker, imageLeft);% before R2016b
[pointsKeyL, isFound] = pointTracker(imageLeft);
% outlier removal
[pointsKeyL, isFound] = outlierRemoval(pointsKeyL,isFound);

pointsNewR = pointsKeyR(isFound,:);
pointsNewL = pointsKeyL(isFound,:);
 
% choose high socres munFeatures for efficiency
if size(pointsNewR, 1)>numFeatures
    randIndex  = randperm(size(pointsNewR, 1),numFeatures);
    pointsKeyR = pointsNewR(randIndex,:); 
    pointsKeyL = pointsNewL(randIndex,:); 
end

% reload tracker
setPoints(pointTracker,pointsKeyR); 
end

function [pointsNewL,pointsNewR,isFoundand] = stereoPointTrack(imageLeft,imageRight,pointTracker) 
% features tacking on new stereo images
% [pointsNewR, isFoundR]  = step(pointTracker, imageRight);   % before R2016b 
[pointsNewR, isFoundR]  = pointTracker(imageRight); 
% [pointsNewL, isFoundL]  = step(pointTracker, imageLeft);% before R2016b
[pointsNewL, isFoundL]  = pointTracker(imageLeft);


% outlier removal
[pointsNewR, isFoundR] = outlierRemoval(pointsNewR,isFoundR);
[pointsNewL, isFoundL] = outlierRemoval(pointsNewL,isFoundL);
% find features in both images
isFoundand = and(isFoundR,isFoundL);
        
end


function S_x = Skew(x)
% Function for Skew symmetric matrix
[row,col] = size(x);    
if ((row == 3)&&(col == 1))
    S_x = [0 -x(3) x(2); x(3) 0 -x(1);-x(2) x(1) 0];
else
    disp('Dimentional error') 
end
end



 

function XhatInv = InvX(Xhat)
Rhat = Xhat(1:3,1:3);
vhat = Xhat(1:3,4);
phat = Xhat(1:3,5); 
XhatInv = [Rhat' -Rhat'*vhat -Rhat'*phat;zeros(2,3) eye(2)];

end

function [points,isFound]=outlierRemoval(points,isFound)
points2 = points(isFound);
pointsmean = mean(points2,1);
temp = points2-pointsmean;
dist = sqrt(diag(temp*temp'));
distmean = mean(dist);
diststd  = std(dist);
index = abs(dist-distmean)<=diststd;
points2= points2(index,:);
isFound = ismember(points(:,1),points2(:,1));
end


function [pointsI,pointsB,pointsNum] = pointsSelection(pointsI,pointsB)
% pointsNum = size(pointsI,2);
% pointsIC   = sum(pointsI,2)/pointsNum; % center of landmarks
% tempI = diag((pointsI - pointsIC)'*(pointsI - pointsIC));
% index = find(tempI<4);
% pointsI = pointsI(:,index);
% pointsB = pointsB(:,index);
% 
% pointsNum = size(pointsB,2);
% pointsBC   = sum(pointsB,2)/pointsNum; % center of landmarks 

% pointsBC  = mean(pointsB,2);
% % pointsBSd = std(pointsB,0,2);
% tempB =  pointsB - pointsBC;
% tempB = sqrt(diag(tempB'*tempB));
% meanB  = mean(tempB);
% stdB = std(tempB);% 
% index = find(abs(tempB-meanB)<max(stdB,0.3));
% pointsI = pointsI(:,index);
% pointsB = pointsB(:,index);

 


 
pointsIC   = mean(pointsI,2); % center of landmarks
pointsBC   = mean(pointsB,2); % center of landmarks
tempI = sqrt(diag((pointsI - pointsIC)'*(pointsI - pointsIC)));
tempB = sqrt(diag((pointsB - pointsBC)'*(pointsB - pointsBC)));
index = abs(tempI-tempB)<=0.01;
% error = abs(tempB - tempI);
% meane = mean(error);
% stde  = std(error);
% index   = find(abs(error-meane)<min(stde,0.009));
pointsI = pointsI(:,index);
pointsB = pointsB(:,index);

pointsNum = size(pointsI,2);
end