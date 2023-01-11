clear all; close all;

syms q0 tht dq0 dtht tau real
syms D L g m phi dt beta ke kd kp k b nPoints s positive
syms d real

q = [q0,tht].'; dq = [dq0,dtht].'; x = [q;dq];
alpha = tht + q0*s;
xsd = D*d*cos(alpha) - L*int(sin(alpha),s,0,s);
ysd = D*d*sin(alpha) + L*int(cos(alpha),s,0,s);

dp = [gradient(xsd,q).';gradient(ysd,q).']; 
M = int(dp.'*dp,d,-.5,.5); 
M = int(M,s,0,1); 
M = m*M;
P = int(xsd*sin(phi) + ysd*cos(phi),d,-.5,.5); 
P = int(P,s,0,1);
P = m*g*P;
G = gradient(P,q);
kappa = [k,0;0,0];
beta = [b,0;0,0];

KE = .5*(dq.'*M*dq);
D = jacobian(KE,dq).'; D = jacobian(D,dq);
syms C real
n=max(size(q));
for k=1:n
    for j=1:n
        C(k,j)=0*g;
        for i=1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i)) + ...
            diff(D(k,i),q(j)) - ...
            diff(D(i,j),q(k)))*dq(i);
        end
    end
end

ddq = M\[tau;0] - M\(C*dq + [b*dq0;0] + [k*q0;0] + G);

A = jacobian([dq;ddq],x);
B = [0;0;0;1];
R = 1;
Q = eye(4);

[D,L,g,m,phi,dt,b,ke,kd,kp,k,nPoints] = getVars();
tht = 0; dtht = 0; dq0 = 0; tau = 0; q0 = 1e-15;
A = subs(A);
A = double(A);
[K,~,~] = lqr(A,B,Q,R);
fprintf("Found LQR Gains To Be:")
disp(K);