function projA=GainMap(A,kR,kv,kp)
A1 = A(1:3,1:3);
% a2 = A(1:3,4);
a3 = A(1:3,5);
PaA1 = (A1-A1')/2;
projA = [kR*PaA1 kv*a3 kp*a3;zeros(2,5)];

end