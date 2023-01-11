syms q0 tht0 q1 tht1 dq0 dtht0 dq1 dtht1 d s real
D = .1; L = 1; %10cm across, 1m long
m = 1; %mass per unit length
phi = 0; %gravity angle
g = 9.8;

x0 = [q0; tht0]; %first member state vector
dx0 = [dq0; dtht0];
x1 = [q1; tht1]; %second memer state vector
dx1 = [dq1; dtht1];
x = [x0; x1];
dx = [dx0; dx1];

%find leg positions
alpha0 = tht0 + s*q0;
alpha1 = tht1 + s*q1; %legs have independent curvatures

xsd0 = D*d*cos(alpha0) - L*int(sin(alpha0),s,0,s);
ysd0 = D*d*sin(alpha0) + L*int(cos(alpha0),s,0,s);

xsd1 = subs(xsd0,s,1) + D*d*cos(alpha1) - L*int((1-s)*sin(alpha1),s,s,0);
ysd1 = subs(ysd0,s,1) + D*d*sin(alpha1) + L*int((1-s)*cos(alpha1),s,s,0);

%get mass matrix
p = [xsd0; ysd0; xsd1; ysd1];
dpdx = jacobian(p,x);
M = int(int( m*dpdx.'*dpdx ,d,-.5,.5),s,0,1);

%get gravity vector
P = int(int( m*g*(xsd0*sin(phi) + ysd0*cos(phi) + xsd1*sin(phi) + ysd1*cos(phi)) ,d,-.5,.5),s,0,1);
G = jacobian(P,x);

%get coriolis matrix
C = sym(zeros(max(size(x))));
n = max(size(x));
for k = 1:n
    for j = 1:n
        for i = 1:n
            C(k,j)=C(k,j)+1/2*(diff(M(k,j),x(i)) + ...
				diff(M(k,i),x(j)) - ...
				diff(M(i,j),x(k)))*dx(i);
        end
    end
end

matlabFunction(M,C,G,"file","genDynamics")
matlabFunction(G,"file","getGroupAction")

T = .5*dx.'*M*dx;
matlabFunction(T,P,"file","getEnergy")

ysd1out = subs(ysd1,[s,d],[1,0]);
matlabFunction(ysd1out,"file","swingFootHeight")