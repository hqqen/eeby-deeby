function [dx] = fdbkCtrlFix(t,x)
if t == 2.252008e+01
    1;
end
q = x(1:2); dq = x(3:4);
[M,C,G,E,PE,KE] = dynMatsFixed(q,dq);
[D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars();

Ed = .5*L*g*m;
% E = .5*(dq.'*M*dq) + .5*(q.'*[k,0;0,0]*q) -(L*g*m*(cos(q(1) + phi + q(2)) - cos(phi + q(2)) + q(1)*sin(phi + q(2))))/q(1)^2;
Ebar = E - Ed;

f = [M\[1;0],...
    -M\(C*dq + [beta*dq(1);0] + [k*q(1);0] + G)];

tau = (-dq(1) + ke*Ebar*beta*dq(1) - kd*f(1,2) - kp*q(1))/(ke*Ebar + kd*f(1,1));

ddq = M\[tau;0] - M\(C*dq + [beta*dq(1);0] + [k*q(1);0] + G);

dx = [dq; ddq];

global tauS; tauS(end+1) = tau;
end