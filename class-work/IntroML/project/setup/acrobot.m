%acrobot dynamics simulator code
%Alex Beyer 2022
%% setup system variables
%setup consistent with Sutton and Barto Ch 11.3 - The Acrobot
syms q1 q2 dq1 dq2 m1 m2 l1 l2 I1 I2 g real %not how inertia works
%get joint and CoM positions
p1 = [l1*sin(q1);...
      -l1*cos(q1)];
pC1 = .5*p1; %CoM position of member 1
p2 = p1 + [l2*sin(q1 + q2);...
      -l2*cos(q1 + q2)];
pC2 = p1 + .5*p2; %CoM position of member 2

%% setup system lagangian
%find kinetic energy
T1 = .5*I1*dq1^2;
T2 = .5*(I2*(dq2^2) + m2*(l1^2)*(dq1^2) + I2*(dq1^2) + 2*m2*l1*(l2/2)*cos(q2)*(dq1^2)) + dq1*dq2*(I2 + m2*l1*(l2/2)*cos(q2));
T = T1 + T2;
%find potential
U = -m1*(l1*.5)*g*cos(q1) - m2*g*(l1*cos(q1) + (.5*l2)*cos(q1 + q2));
L = T - U; %lagrangian
H = T + U; %hamiltonian
%% get equations of motion in standard form
q = [q1,q2].'; dq = [dq1,dq2].';
M = jacobian(jacobian(T,dq).',dq); %inertia matrix
G = jacobian(U,q); %gravity vector
% Cdq = jacobian(M*dq,q)*dq -.5*(jacobian(M*dq,q)*dq).'*dq
% n = max(size(q));
% C = sym(zeros(n,n));
% for k = 1:n
%     for j = 1:n
%         for i = 1:n
%             C(k,j)=C(k,j)+1/2*(diff(M(k,j),q(i)) + ...
% 				diff(M(k,i),q(j)) - ...
% 				diff(M(i,j),q(k)))*dq(i);
%         end
%     end
% end
% C = simplify(C);

%The way in which the textbook (and presumably the gym model) find the
%Coriolis matrix is nonstandard; to ensure my dynamics stay consistent I'm
%going to hardcode it
C = [-2*m2*l1*(l2/2)*sin(q2)*dq2, -m2*l1*(l2/2)*sin(q2)*dq2;...
    m2*l1*(l2/2)*sin(q2)*dq1,    0]; %coriolis matrix

B = [0;1]; %control mapping

% matlabFunction(M,C,G,'file', 'acrobotDyn');
% matlabFunction(H,L,T,U,'file','systemEnergy');
% matlabFunction(p1,p2,'file','getJointPosition');
