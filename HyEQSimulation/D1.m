function [inside] = Dset(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                
%
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A1 A2
%% state
% xIn = [x0;xhat0;vP10;xhat0;vP20;tau;time];  % 15+(15+3)+(15+45) +2
% % real system states
% vR     = x(1:9,1);
% p      = x(10:12,1);
% v      = x(13:15,1);
% % A-HINO states 
% index  = 15;
% lP1    = length(A1)*(length(A1)+1)/2;
% vRhat1 = x(index+1:index+9,1);
% phat1  = x(index+10:index+12,1);
% vhat1  = x(index+13:index+15,1);
% vP1    = x(index+16:index+15+lP1,1);
% % iEKF states
% index  = index+15+lP1;
% lP2    = length(A2)*(length(A2)+1)/2;
% vRhat2 = x(index+1:index+9,1);
% phat2  = x(index+10:index+12,1);
% vhat2  = x(index+13:index+15,1);
% vP2    = x(index+16:index+15+lP2,1);
% % time
tau    = x(end-1,1);
% time   = x(end,1); 

if (tau<=0)
    inside = 1;
else
    inside = 0;
end


end