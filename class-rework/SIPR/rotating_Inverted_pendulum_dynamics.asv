clc

syms a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20
syms da0 da1 da2 da3 da4 da5 da6 da7 da8 da9 da10 da11 da12 da13 da14 da15 da16 da17 da18 da19 da20
syms s u v d D L phi m g beta k tau1 theta dtheta

a_vect = [a0;a1;a2;a3;a4;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14;a15;a16;a17;a18;a19;a20];
da_vect = [da0;da1;da2;da3;da4;da5;da6;da7;da8;da9;da10;da11;da12;da13;da14;da15;da16;da17;da18;da19;da20];
v_vect = [1;v;v^2;v^3;v^4;v^5;v^6;v^7;v^8;v^9;v^10;v^11;v^12;v^13;v^14;v^15;v^16;v^17;v^18;v^19;v^20];

% phi = 0;
% D = 0.1;
% L = 1;
% m = 1;
% g =9.81;
% syms k
% beta = 0.1;
% syms tau1

i = 1; %% order

aa = a_vect(1:i);
vv = v_vect(1:i);
daa = da_vect(1:i);

q = aa.'*vv;

alpha_ = int(q,v,0,v) + theta;

xx = sin(alpha_);

x_ = -L*int(xx,v,0,v);

yy = cos(alpha_);

y_ = L*int(yy,v,0,v);

x = x_ + d*D*cos(alpha_);
y = y_ + d*D*sin(alpha_);

% RR = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% new_xy = RR*[x;y];

% x = new_xy(1);
% y = new_xy(2);

vars = [aa;theta];
dvars = [daa;dtheta];
jsd = (jacobian([x;y],vars));

% rho = m*dirac(v-1);
rho = m;

B1 = rho*(jsd.')*jsd;

B_mtx1 = int(B1,d,-0.5,0.5);

D_mtx = int(B_mtx1,v,0,1);
D_mtx1= simplify(expand(D_mtx))

% 
% %%Kinetic Energy%%
vx = diff(x,a0)*da0+ diff(x,theta)*dtheta;
vy = diff(y,a0)*da0+ diff(y,theta)*dtheta;
KE1 = int(0.5*rho*(vx^2 + vy^2),d,-0.5,0.5);
KE = int(KE1,v,0,1);
%%Potential Energy%%

PE = int(int((rho*g*(sin(phi)*x + cos(phi)*y)),d,-0.5,0.5),v,0,1);


Total_Energy = KE + PE + int(0.5*k*q*q,v,0,1);
%%Solving Lagrangean%%

D_mtx2=(jacobian(KE,[daa;dtheta]).');
D_mtx2=(jacobian(D_mtx2,[daa;dtheta]));
D_mtx2 = simplify(expand(D_mtx2))

syms C_mtx real
n=max(size(vars));
for k=1:n
	for j=1:n
		C_mtx(k,j)=0*g;
		for i=1:n
			C_mtx(k,j)=C_mtx(k,j)+1/2*(diff(D_mtx(k,j),vars(i)) + ...
				diff(D_mtx(k,i),vars(j)) - ...
				diff(D_mtx(i,j),vars(k)))*dvars(i);
		end
	end
end
C_mtx=simplify(C_mtx)


syms k u
% g_field = rho*g*sin(phi)*x + rho*g*cos(phi)*y;
% 
% G1 = int(g_field,d,-0.5,0.5);
% 
% G = int(G1,v,0,1);
% 
% G_vect = jacobian(G,aa).'

G_vect = jacobian(PE,vars).'

H = 1;

Tau=tau1;


K = [k 0; 0 0];
Damp = [beta 0; 0 0];


hqq  = C_mtx*[daa;dtheta] + G_vect + K*[aa;theta] + Damp*[daa;dtheta];

B_tilde = D_mtx;


% ddqq = inv(B_tilde)*(Tau - hqq);


% new_tau = (hqq(1) - (B_tilde(2,1)/B_tilde(2,2))/hqq(2)) + (B_tilde(1,1) - B_tilde(2,1)*B_tilde(2,1)/B_tilde(2,2))*u;

xdoubledot=[da0; dtheta; inv(D_mtx)*([1;0]*Tau - C_mtx*[daa;dtheta] - G_vect - K*[aa;theta] -Damp*[daa;dtheta])];
% xdoubledot=[da0; dtheta; inv(D_mtx)*([0;1]*Tau - C_mtx*[daa;dtheta] - G_vect - K*[aa;theta] -Damp*[daa;dtheta])];
fx=[xdoubledot];
XX=[aa;theta;da0; dtheta];
gx=[XX];
% gx=simplify([inv(Jacob.')*(-B*ddtheta-(C+Damp)*dtheta-G-K*theta+Tau)]);

%feedback linearized
% fx = [da0;dtheta;tau1;inv(B_tilde(2,2))*(-hqq(1) - B_tilde(1,2)*tau1)];
% 
% fx = [da0;dtheta;inv(B_tilde(2,1))*(-hqq(2) - B_tilde(2,2)*tau1);tau1];
% Tau=u;

% 
A_matrix=jacobian(fx,XX)
B_matrix=jacobian(fx,Tau)
C_matrix=jacobian(gx,XX);
D_matrix=jacobian(gx,Tau);
 %%

ref_theta=[1e-15,1e-15];
ref_dtheta=[1e-16,1e-16];
ref_ddtheta=[0,0];
mass=2;
Len=1;
k_ref=0.5;
beta_ref=0.1;
ref_tau=0;
D_ref = 0.1;
g_ref=9.81;
phi_ref =0;

ref_point=[ref_theta,ref_dtheta,mass,Len,D_ref,k_ref,beta_ref,ref_tau,g_ref,phi_ref];

AA=double(subs(A_matrix,[a0,theta,da0,dtheta,m,L,D,k,beta,tau1,g,phi],ref_point));
BB=double(subs(B_matrix,[a0,theta,da0,dtheta,m,L,D,k,beta,tau1,g,phi],ref_point));
CC=double(subs(C_matrix,[a0,theta,da0,dtheta,m,L,D,k,beta,tau1,g,phi],ref_point));
DD=double(subs(D_matrix,[a0,theta,da0,dtheta,m,L,D,k,beta,tau1,g,phi],ref_point));


System=ss(AA,BB,CC,DD)

% Checking controllability
con_rank = rank(ctrb(AA,BB))

if con_rank==4
    disp('Controllable.... Hurray!!');
else
    disp('UNCONTROLLABLE !!');
end

% Checking observability
ob_rank = rank(obsv(AA,CC))

if ob_rank==4
    disp('Observable.... Hurray!!');
else
    disp('UNOBSERVABLE !!');
end


Q = [1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1];

lqr(System,Q,1)
