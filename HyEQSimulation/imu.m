function [a,omega]=imu(R,t)
global w g e3
omega = 1*[sin(0.3*pi*t) 0 1]';
a     = R'*(-10*w^2*[cos(w*t) sin(w*t) 0]'+g*e3);

end