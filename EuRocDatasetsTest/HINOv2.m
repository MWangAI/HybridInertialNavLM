clear all
close all
clc

% Navigation on SE_2(3) with biased IMU
% Limited number of featrues in the pointTracker!
% for example 30. This method saves time for the feature tracking
%% Read dataset

% Load stereo camera intrinsic parameters (calibrated using MATLAB)
load('stereoParamsdataset.mat'); 

% Set the path of different experiments
setpath = 1;
switch  setpath
    case 1
        path        = 'E:\Dataset\V1_01_easy'; 
    case 2
        path        = 'E:\Dataset\V1_02_medium'; 
    case 3
        path        = 'E:\Dataset\V2_01_easy'; 
    case 4
        path        = 'E:\Dataset\V2_02_medium'; 
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
numFeatures = 100;
index       = 1;
pointTracker = vision.PointTracker('BlockSize',[11 11],'MaxBidirectionalError',1); %Reduce the blocksize can improve the accuracy
[imageLeft,imageRight]                       = stereoRead(path,CamRead,index,stereoParams);
[pointsKeyR,pointsKeyL,pointTracker,isFound] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker);
% Show matched features
% figure; ax = axes;
% showMatchedFeatures(imageLeft,imageRight,pointsKeyL(isFound,:),pointsKeyR(isFound,:),'montage','Parent',ax);
% % title(ax, 'Candidate point matches');
% legend(ax, 'Matched points left','Matched points right');
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','Features.eps')
%% Initialize observer

g  = 9.81;
e3 = [0,0,1]';
 
P  = diag([1*ones(1,3) 1*ones(1,3) 0.0001*ones(1,3)]); 
kw = 30*0.0050;
Q  = 0.005*eye(3);
V  = 0.005*diag([1*ones(1,3) 1*ones(1,3) 1*ones(1,3)]); 
    
    
p  = Groundtruth(1,2:4)';
R  = quat2rotm(Groundtruth(1,5:8));
v  = Groundtruth(1,9:11)';
 



Rkey = R;
pkey = p;
vkey = v; 
points3D  = funPoints3D(pointsKeyL,pointsKeyR,isFound,stereoParams,TB2C0);
pointsKey = zeros(3,length(isFound)); 
pointsKey(:,isFound==1)   = Rkey*points3D+pkey;  %inertial frame landmarks  

u      = e3;%randn(3,1);
theta  = 0.99*pi;%(0.5-rand(1))*0.4*pi;
Rq     = expm(Skew(u)*theta);
Xhat   = [Rq*R 0*v+0.0*randn(3,1) 0*p+0.0*randn(3,1);
         zeros(2,3) eye(2)];
% Xhat   = [Rq*R 1*v+0.0*randn(3,1) 1*p+0.0*randn(3,1);
%          zeros(2,3) eye(2)];
% Xhat   = [eye(3) 0*v 0*p;zeros(2,3) eye(2)];

bwhat  = zeros(3,1);
bahat  = zeros(3,1);
 


%% Loop
Tfinal     = ceil(length(TimerIMU)*0.95);%20000;%
PGrdtruth = Groundtruth(1:Tfinal,2:4)';
Pestimate  = zeros(3,Tfinal); 
normR      = zeros(1,Tfinal);
normV      = zeros(1,Tfinal); 
normbw     = zeros(1,Tfinal);
normba     = zeros(1,Tfinal); 
keytime        = zeros(size(CamRead));
keytime(index) = 1;
distemp = 0;
for k=1:1:Tfinal
    if (mod(floor(100*k/Tfinal),10)==0)&&(floor(100*k/Tfinal)>distemp)
    clc
    disp(['Complete  ' num2str(floor(100*k/Tfinal)) '%'])
    distemp = floor(100*k/Tfinal);
    end
    % errors
    p       = Groundtruth(k,2:4)';
    R       = quat2rotm(Groundtruth(k,5:8));
    v       = Groundtruth(k,9:11)';     
    normR(1,k)      = trace(eye(3)-R*Xhat(1:3,1:3)')/4;  
    normV(1,k)      = norm(v-Xhat(1:3,4));  
    normbw(1,k)     = norm(bw'-bwhat);
    normba(1,k)     = norm(ba'-bahat);
    Pestimate(:,k)  = Xhat(1:3,5)';  
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
    vhat = vhat + (-g*e3+Rhat*(a-bahat))*dT;
    Xhat = [Rhat vhat phat;zeros(2,3) eye(2)];
 
    
 
    
    A = [-Skew(omega-bwhat) eye(3) zeros(3);zeros(3) -Skew(omega-bwhat) eye(3);zeros(3,9)];
    C = [eye(3) zeros(3) zeros(3)];
    
%     approach 1
%     Ad    = expm(A*dT);
%     P     = Ad*P*Ad'+ 0.1*V; %
    
%     apprach 2
    P = RK4(A,P,V,dT);% 
%     or
%     P  = P + (A*P + P*A'+V)*dT;
% %     
%     HH = [-A' zeros(6);V A];
%     MM = (expm(HH*dT)-eye(12))/dT;
%     dP = (MM(7:12,7:12)*P - P*MM(1:6,1:6)+ MM(7:12,1:6))/(eye(6)+dT*MM(1:6,1:6));
%     P3 = P + dP*dT;
     
    
    
    if TimerIMU(k+1)>= TimerCam(index+1,1)
        index = index +1;               
        [imageLeft,imageRight]             = stereoRead(path,CamRead,index,stereoParams);
        [pointsNewL,pointsNewR,isFoundand] = stereoPointTrack(imageLeft,imageRight,isFound,pointTracker,pointsKeyR);         
        
%         savefigures(imageLeft,imageRight,pointsNewL,pointsNewR,isFoundand,PGrdtruth,Pestimate,index-1,k)
        
%         if (mod(index-2,2)==0)
%             ii = floor((index-2)/2)+1;
%             savefigures(imageLeft,imageRight,pointsNewL,pointsNewR,isFoundand,PGrdtruth,Pestimate,ii,k)
%         end
        
        
        pointsI   = pointsKey(:,isFoundand==1);
%         pointsNewL = pointsNewL + 5*randn(size(pointsNewL)); %add noise for analysis
%         pointsNewR = pointsNewR + 5*randn(size(pointsNewR));
        pointsB   = funPoints3D(pointsNewL,pointsNewR,isFoundand,stereoParams,TB2C0);  
        [pointsI,pointsB] = pointsSelection(pointsI,pointsB); 
        pointsNum = size(pointsI,2);
    

        
        if (pointsNum>=3)
            % 
            Kn      = eye(pointsNum)/pointsNum;
            pointsC = sum(pointsI*Kn,2); % center of landmarks
            M       = (pointsI-pointsC)*Kn*(pointsI-pointsC)';
            Mbar   = (trace(M)*eye(3)-M)/2;
            if (max(eig(Mbar))>0.01) 
                kR    = 0.4/max(eig(Mbar));
            else
                kR    = 20;
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
            Kp     = Rhat*KP(1:3,1:3)*Rhat'/sum(diag(Kn));%Rhat*KP(1:3,1:3)*Rhat';%
            Kv     = Rhat*KP(4:6,1:3)*Rhat'/sum(diag(Kn)); %Rhat*KP(4:6,1:3)*Rhat';%   
            Ka     = Rhat*KP(7:9,1:3)*Rhat'/sum(diag(Kn)); %Rhat*KP(7:9,1:3)*Rhat';%   
            Temp1    = XcInv*(r-Xhat*b)*Kn*r'*XcInv'; 
            Temp2    = GainMap(Temp1,kR,Kv,Kp); 
            Delta   = Xc*Temp2*XcInv;
 
           
%             Rhat    = expm(Delta(1:3,1:3))*Rhat;
%             phat    = phat + Delta(1:3,1:3)*(phat-pointsC) + Delta(1:3,5);
%             vhat    = vhat + Delta(1:3,1:3)*vhat + Delta(1:3,4);
%             Xhat    = [Rhat vhat phat;zeros(2,3) eye(2)];
            XhatInv = InvX(Xhat); 
            Xhat    = Xhat*expm(XhatInv*Delta*Xhat); % matrix form
            psiR   =  [Delta(3,2) Delta(1,3) Delta(2,1)]';
            bwhat = bwhat +  (-kw*Rhat'*psiR);
            bahat = bahat + (-Rhat'*Ka*Temp1(1:3,5));
            
            P       = (eye(size(P))-KP*C)*P; 
            
        end %endif
        % Creat new Keyframe if pointsNum is less than 6     
        minFeatures = 6;
        if pointsNum<=minFeatures %size(find(isFoundand==1), 1) < minFeatures       
            [pointsKeyR,pointsKeyL,pointTracker,isFound] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker);                 
            pkey     = Groundtruth(k,2:4)';
            Rkey     = quat2rotm(Groundtruth(k,5:8));
            points3D = funPoints3D(pointsKeyL,pointsKeyR,isFound,stereoParams,TB2C0);
            pointsKey               = zeros(3,length(isFound)); 
            pointsKey(:,isFound==1) = Rkey*points3D+pkey;  %inertial frame landmarks  
            keytime(index) = 1;

         end %endif
    end %endif
    
    
    
end %endfor
 
 
release(pointTracker);
 
clearvars -except TimerIMU PGrdtruth Tfinal Pestimate TimerCam keytime normR normbw normV normba
%% Plot

errorP     = Pestimate - PGrdtruth;
errorPnorm = sqrt(sum(errorP.^2,1));
RMSE       = sqrt(mean(errorPnorm.^2)) 

 

Time = (TimerIMU(1:Tfinal)-TimerIMU(1));


figure 
subplot(2,2,1)
plot(Time,sqrt(normR(1,:))','b-','linewidth',1), grid on
ylabel('$|\tilde{R}|_I$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
xlim([0 Time(end)])
ylim([0 1.0])
    axes('position',[.22 .62 .20 .22])
    box on % put box around new pair of axes
%     indexOfInterest1  = (Time<=1); 
    indexOfInterest1  = (Time>=0.95*Time(end)) & (Time<=Time(end)); 
    plot(Time(indexOfInterest1),sqrt(normR(1,indexOfInterest1))','b-','linewidth',1), grid on
    axis tight
    
subplot(2,2,2) 
plot(Time,errorPnorm,'b-','linewidth',1), grid on
xlim([0 Time(end)])
% ylim([0 0.9]) 
ylabel('$\|p-\hat{p}\|$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
    axes('position',[.72 .62 .20 .22])
    box on % put box around new pair of axes
    indexOfInterest1  = (Time<=5); 
%     indexOfInterest1  = (Time>=0.98*Time(end)) & (Time<=Time(end)); 
    plot(Time(indexOfInterest1), (errorPnorm(indexOfInterest1))','b-','linewidth',1), grid on
    axis tight
 
subplot(2,2,3)
plot(Time,normV(1,:),'b-','linewidth',1), grid on
ylabel('$\|v-\hat{v}\|$','interpreter','latex')
xlabel('$t(s)$','interpreter','latex')
xlim([0 Time(end)])
    axes('position',[.22 .22 .20 .22])
    box on % put box around new pair of axes
    indexOfInterest1  = (Time<=5); 
%     indexOfInterest1  = (Time>=0.98*Time(end)) & (Time<=Time(end)); 
    plot(Time(indexOfInterest1), (normV(1,indexOfInterest1))','b-','linewidth',1), grid on
    axis tight
% ylim([0 1.6])  

subplot(2,2,4)
yyaxis left
plot(Time,normbw(1,:),'b-','linewidth',1), grid on
ylabel('$\|b_\omega-\hat{b}_\omega\|$','interpreter','latex')
yyaxis right
plot(Time,normba(1,:),'r-','linewidth',1), grid on
ylabel('$\|b_a-\hat{b}_a\|$','interpreter','latex')
xlim([0 Time(end)])
xlabel('$t(s)$','interpreter','latex')

% ylim([0 0.16]) 
% legend('HINO', 'HINO-CRE2')
%  set(gcf, 'Renderer', 'Painters');
%  print('-depsc','SimulationError3.eps')
 
 
% figure
% plot(TimerCam-TimerIMU(1),Numlandmark,'.'), hold on
% xlim([0 Time(end)])
% Keyframetime = TimerCam(keytime==1);  
% plot(Keyframetime-TimerIMU(1),zeros(size(Keyframetime)),'*')



figure 
plot3(PGrdtruth(1,:), PGrdtruth(2,:),PGrdtruth(3,:),'g--','linewidth',1), hold on
plot3(Pestimate(1,:), Pestimate(2,:),Pestimate(3,:),'b-','linewidth',1) 
% title(['RMSE=' num2str(RMSE)])
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
grid on
legend('True trajectory','HINO-CRE')
 set(gcf, 'Renderer', 'Painters');
 print('-depsc','3dTrajectory3.eps')
 

 function savefigures(imageLeft,imageRight,pointsKeyL,pointsKeyR,isFound,PGrdtruth,Pestimate,ii,k)
%     figure
    h1=subplot(2,1,1); 
    h1.Position = h1.Position + [-0.1 -0.05 0.15 0.1];
    [~,col] = size(imageLeft);
    I = [imageLeft imageRight];
%     I = insertMarker(I,pointsKeyL(isFound,:),'o');
    imshow(I), hold on;
    plot(pointsKeyL(isFound,1),pointsKeyL(isFound,2),'yo')
    plot(col+pointsKeyR(isFound,1),pointsKeyL(isFound,2),'y+')
    
%     showMatchedFeatures(imageLeft,imageRight,pointsKeyL(isFound,:),pointsKeyR(isFound,:),'montage','Parent',ax);
    h2=subplot(2,1,2);
    h2.Position =  [0.20 0.08 0.58 0.48];
    plot3(PGrdtruth(1,:), PGrdtruth(2,:),PGrdtruth(3,:),'g--','linewidth',1), hold on
    plot3(Pestimate(1,1:k), Pestimate(2,1:k),Pestimate(3,1:k),'b-','linewidth',1.5) 
    % title(['RMSE=' num2str(RMSE)])
%     axis([-2 2 -2.5 4 -2 2])
    minp = min(PGrdtruth')-0.3;
    maxp = max(PGrdtruth')+ 0.3;
    axis([minp(1) maxp(1)  minp(2) maxp(2) minp(3) maxp(3)])
    legend('True trajectory','HINO-CRE2','Position',[0.67 0.46 0.22 0.08]);
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    grid on    
    view(104.69, 66.61)
    
    
    saveas(gcf,['Video2\' num2str(ii,'%04.f') '.png'])
 end
 
 

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

function [positionsB] = funPoints3D(pointsLeftFrame,pointsRightFrame,isFoundand,stereoParams,TB2C0)
pointsFeaturesR  = pointsRightFrame(isFoundand, :);
pointsFeaturesL  = pointsLeftFrame(isFoundand, :);
%     position on the frame of left cam 
positionCam0 = 0.001*triangulate(pointsFeaturesL,pointsFeaturesR,stereoParams);
positionCam0 = positionCam0';
%     positionCam0 = R'(positionsB-p)
R = TB2C0(1:3,1:3);
p = TB2C0(1:3,4);
positionsB = R*positionCam0 + p;    %R'*(positionCam0-p);%  
    
end

function [pointsKeyR,pointsKeyL,pointTracker,isFound] = funPointTracker(imageLeft,imageRight,numFeatures,pointTracker)

% feature detection on the right image of keyframe and set the points as
% pointTracker
 

points          = detectMinEigenFeatures(imageRight);
pointsKeyR      = points.Location;

% randomly choose
if size(pointsKeyR, 1)>numFeatures
    randIndex  = randperm(size(pointsKeyR, 1),numFeatures);
    pointsKeyR = pointsKeyR(randIndex,:); 
end

release(pointTracker);
% Initialize the tracker with the initial point locations and the initial video frame.
initialize(pointTracker, pointsKeyR, imageRight);
% track features on the left image of keyframe
% [pointsKeyL, isFound] = step(pointTracker, imageKeyframeL);% before R2016b
[pointsKeyL, isFound] = pointTracker(imageLeft);
% outlier removal
[isFound] = outlierRemoval(pointsKeyL,pointsKeyR,isFound,50,2);
end

function [pointsNewL,pointsNewR,isFoundand] = stereoPointTrack(imageLeft,imageRight,isFound,pointTracker,pointsKeyR) 
% features tacking on new stereo images
% [pointsNewR, isFoundR]  = step(pointTracker, imageRight);   % before R2016b 
[pointsNewR, isFoundR]  = pointTracker(imageRight); 
% [pointsNewL, isFoundL]  = step(pointTracker, imageLeft);% before R2016b
[pointsNewL, isFoundL]  = pointTracker(imageLeft);


% outlier removal
[isFoundR] = outlierRemoval(pointsNewR,pointsKeyR,isFoundR,60,2);
[isFoundL] = outlierRemoval(pointsNewL,pointsKeyR,isFoundL,50,2);
% find features in both images
isFoundand = and(isFound,and(isFoundR,isFoundL));
        
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

function [isFound]=outlierRemoval(points1,points2,isFound,S,D)
for i=1:2
    points = points1(isFound,i)-points2(isFound,i);
    points1New = points1(isFound,i);
    index  = abs(points)<=S;
    points = points(index);
    points1New = points1New(index);

    distmean = mean(points);
    diststd  = std(points);
    index = abs(points-distmean)<=max(diststd,D);
    points1New = points1New(index);
    isFound = ismember(points1(:,i),points1New); 
end
end


function [pointsI,pointsB] = pointsSelection(pointsI,pointsB)


 
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
 
end