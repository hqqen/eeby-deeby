clear all; close all;
%symbolically calulate dynamic matrices M (inertia), C (coriolis) and
%G(gravity) then find linearized system matrices A & B for LQR upright
%stabilization problem

%mode 1 - keep everything symbolic for as long as possible
%mode 2 - substitute values in as early as possible (debug to test for
%more matlabFunction bugs)
mode = 1;
%nSimplify is the number of simplification passes to do on the matrices,
%try a range of them to see what works best
nSimplify = 1;
%linearizePt is the value of q0 to linearize about since doing it at 0
%causes problems
linearizePt = 1e-30;


%setup vars
% if mode == 1
    syms q0 tht s D d L l m phi g t ye beta dq0 dtht tau real
    assumeAlso([s,D,L,m,g,t,beta],'positive');
% %     assumeAlso(tht > 0 & tht < 2*pi);
% elseif mode == 2
%     [D,L,g,m,phi,dt,beta,~,~,~,~,nPoints] = getVars();
%     syms q0 tht dq0 dtht d s l t ye tau real
%     assumeAlso([s],'positive');
%     assumeAlso(tht > 0 & tht < 2*pi);
% end
q = [q0; tht]; dq = [dq0; dtht]; x = [q; dq];

%% Get dynamic matrices
%get position along member
alpha = tht + int(q0,l,0,s);
xsd = D*d*cos(alpha) - L*int(sin(alpha),s,s,0);%<-- the order of these variable matters???????
ysd = D*d*sin(alpha) + L*int(cos(alpha),s,s,0);
dp = [gradient(xsd,q),gradient(ysd,q)].'; %<---- CHECK TRANSPOSITIONS 

%get inertia matrix M
M = int(int(m*(dp.'*dp),d,-.5,.5),s,0,1);
% M = simplify(M,nSimplify);
fprintf("Inertia Matrix Found!\n")

%find potential energy
P = int(int(m*g*(xsd*sin(phi) + ysd*cos(phi)),d,-.5,.5),s,0,1);
% P = simplify(P,nSimplify);

%get gravity vector
G = gradient(P,q);
% G = simplify(G,nSimplify);
fprintf("Gravity Vector Found!\n")

%get coriolis matrix via christoffel symbols
C = sym(zeros(size(M)));
K = .5*(dq.'*M*dq); D = jacobian(jacobian(K,q).',q);
for k = 1:size(q,1)
    for j = 1:size(q,1)
        C(k,j) = 0*g;
        for i = 1:size(q,1)
            C(k,j) = C(k,j) + .5*(diff(M(j,k),q(j)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*dq(i);
%             C(k,j) = C(k,j) + .5*(diff(D(j,k),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)))*dq(i);
            fprintf("Coriolis Matrix Coordinate (" + k + "," + j + ") Iterant " + i + " Done\n")
        end
    end
end
fprintf("Coriolis Matrix Done!\n")
% Cdq = jacobian(M*dq,q)*dq - .5*(jacobian(M*dq,q).')*dq;
% C = simplify(C,nSimplify);
% C = jacobian(M*dq,q)*dq - .5*jacobian(M*dq,q).';

%find MPFL normal form
% ddqTau = M\[1;0];
% ddq =  -M\(C*dq + [beta*dq(1);0] + [k*q(1);0] + G);
% f = [ddqTau(1), ddq(1);...
%      ddqTau(2), ddq(2)];

%% Find system mats for LQR stabilization
%if my LQR gains are off that means the issue is in the dynamics
%get system matrix A
%dx = Ax + Bu
dx = [...
    dq;
    M\[tau;0] - M\(C*dq + [beta*dq0;0] + [k*q0;0] + G)];
A = jacobian(dx, x); % this is wrong, the derivative for ddq w/r.t. dq0 and dtht probably shouldnt be 0
fprintf("Linearized System Matrix Found!\n")
% if mode == 1
    [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars();
    q0 = linearizePt; tht = 0; dq0 = 0; dtht = 0; tau = 0;
    A = double(subs(A));
% elseif mode == 2
%     q0 = linearizePt; tht = 0; dq0 = 0; dtht = 0; tau = 0;
%     A = double(subs(A));
% end
%define control matrix B - should only map inputs to ddq0 (dx(3))
B = [0;0;0;1];
R = 1;
Q = eye(4);
[K,~,~] = lqr(A,B,Q,R);
fprintf("Found LQR Gains To Be:")
disp(K);
%% convert to functions
% matlabFunction(M,C,G,P,"file","dynMats");%,"Optimize",false)
% matlabFunction(x,y,"file","curvatureMap")
% matlabFunction(A,B,"file","sysMatsNoSub","Optimize",false)
