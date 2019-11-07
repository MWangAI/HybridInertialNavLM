
function P=RK4(A,P,V,dT)
DP1 = (A*P+P*A'+V);
P1  = P + DP1*dT/2;
DP2 = (A*P1+P1*A'+V);
P2  = P + DP2*dT/2;
DP3 = (A*P2+P2*A'+V);
P3  = P + DP3*dT;
DP4 = (A*P3+P3*A'+V);
P   = P + dT*(DP1+2*DP2+2*DP3+DP4)/6;
end
