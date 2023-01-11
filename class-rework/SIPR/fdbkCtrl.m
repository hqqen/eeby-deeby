function [dx] = fdbkCtrl(t,x)
%integrator to pass into ode45
%takes time t and state vector x = [q0,tht,dq0,dtht]

%setup vars
q = x(1:2); %q(2) = sign(q(2))*mod(q(2),pi); %tht in S1
% if q(2) > 2*pi
%     q(2) = q(2) - 2*pi;
% elseif q(2) < 0
%     q(2) = q(2) + 2*pi;
% end
% if q(1) == 0
%     q(1) = 1e-30;
% end
dq = x(3:4);
[M,C,G,P] = dynMats(q,dq);
[~,L,g,m,~,~,beta,ke,kd,kp,k,~] = getVars();

%find control input
%rewrite control equation in MPFL normal form
%find MPFL normal form
ddqTau = M\[1;0];
ddq =  -M\(C*dq + [beta*dq(1);0] + [k*q(1);0] + G);
f = [ddqTau(1), ddq(1);...
     ddqTau(2), ddq(2)];
% E = .5*dq.'*M*dq + .5*k*q(1)^2 + ((L*g*m)/(q(1)^2))*(cos(q(2)) - q(1)*sin(q(2)) - cos(q(1) - q(2)));
E = .5*(dq.'*M*dq) + .5*(q.'*[k,0;0,0]*q) + P;
Ed = .5*m*g*L;
Estar = E - Ed;
tau = (-dq(1) + ke*Estar*beta*dq(1) - kd*f(1,2) - kp*q(1))/(ke*Estar + kd*f(1,1));
% tau = 0;
global tauS; tauS(end+1) = tau;
%calculate state upadte equation
% clear ddq;
% ddq = M*[tau;0] - M\((C*dq) + ([beta*dq(1);0]) + ([k*q(1);0]) + G);

%make new state update vector
dx = [...
     dq;
     f(1,1)*tau + f(1,2);
     f(2,1)*tau + f(2,2)];
end