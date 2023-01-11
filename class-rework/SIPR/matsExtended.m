close all; clear all;
syms q0 dq0 q1 dq1 tht0 dtht0 tht1 dtht1 s m L d D k beta tau0 tau1 b g phi footX footY dfootX dfootY hipX hipY real

%leg 0 - stance
%leg 1 - swing

q = [q0 tht0 q1 tht1].'; dq = [dq0 dtht0 dq1 dtht1].'; x = [q;dq];
qe = [q; footX; footY]; dqe = [dq; dfootX; dfootY]; xe = [qe;dqe];
tau = [tau0; tau1];

%define position along leg; note that swing leg (leg 1) posn is a fn of
%stnace leg (leg 0) posn
%also find first time derivative of the position via chain rule, then get
%the jacobian of leg positions w/r.t. state vector q

alpha0 = int(q0,s,0,s) + tht0;
xsd0 = d*D*cos(alpha0) - L*int(sin(alpha0),s,0,s);
ysd0 = d*D*sin(alpha0) + L*int(cos(alpha0),s,0,s);
dxsd0 = diff(xsd0,q0)*dq0 + diff(xsd0,tht0)*dtht0;
dysd0 = diff(ysd0,q0)*dq0 + diff(ysd0,tht0)*dtht0;

alpha1 = int(q1,s,0,s) + tht1; %dont add alpha0 to this; the legs are connected @ hip by ideal joint
xsd1 = d*D*cos(alpha1) - L*int(sin(alpha1),s,0,s) + hipX;% + subs(xsd0,s,1); %leave this as just the foot position
ysd1 = d*D*sin(alpha1) + L*int(cos(alpha1),s,0,s) + hipY;% + subs(ysd0,s,1);
dxsd1 = diff(xsd1,q1)*dq1 + diff(xsd1,tht1)*dtht1;
dysd1 = diff(ysd1,q1)*dq1 + diff(ysd1,tht1)*dtht1;

psd = jacobian([xsd0;ysd0;xsd1;ysd1],q);

%find inertia matrix M
M = m*(psd.'*psd);
M = int(M,d,-.5,.5);
M = int(M,s,0,1);
M = simplify(expand(M));

PE = int(int((m*g*(sin(phi)*xsd0 + cos(phi)*ysd0 + sin(phi)*xsd1 + cos(phi)*ysd1)),d,-0.5,0.5),s,0,1);
G = jacobian(PE,q).';

n = max(size(q));
for k = 1:n
    for j = 1:n
        C(k,j) = 0*g;
        for i = 1:n
            C(k,j)=C(k,j)+1/2*(diff(M(k,j),q(i)) + ...
				diff(M(k,i),q(j)) - ...
				diff(M(i,j),q(k)))*dq(i);
        end
    end
end
C = simplify(C);
syms k 
kappa = [k 0 0 0; 0 0 0 0; 0 0 k 0; 0 0 0 0];
b = [beta 0 0 0; 0 0 0 0; 0 0 beta 0; 0 0 0 0];

ddx = [dq;...
       M\([tau0;0;tau1;0] - C*dq - b*dq - kappa*q - G)];
A2 = jacobian(ddx,x);
B2 = jacobian(ddx,tau);
C2 = jacobian(x,x);
D2 = jacobian(x,tau);

[D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars();
dq0 = 0; tht0 = 1e-15; dtht0 = 1e-15; tau0 = 0; q0 = 1e-30;
dq1 = 0; tht1 = 1e-15; dtht1 = 1e-15; tau1 = 0; q1 = 1e-30;
A2 = double(subs(A2));
B2 = double(subs(B2));
C2 = double(subs(C2));
D2 = double(subs(D2));

sys = ss(A2,B2,C2,D2);
Q = eye(8);

lqr(sys,Q,1)