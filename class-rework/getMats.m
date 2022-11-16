close all; clear all;
syms q0 dq0 q1 dq1 tht0 dtht0 tht1 dtht1 s m L d D k beta tau0 tau1 b g phi footX footY dfootX dfootY hipX hipY dhipX dhipY real

%script to symbolically calculate matrices needed for soft bipedal gait
%generation & convert to function
%two soft members pinned at the "hip"; stance foot pos'n is stored btwn
%foot strikes with both legs having base angle msrd from hip
%q0/q1 being respective curvatures
%leg 0 - stance
%leg 1 - swing

%def coordinate vectors
q = [q0 tht0 q1 tht1].'; dq = [dq0 dtht0 dq1 dtht1].'; x = [q;dq];
qe = [q; footX; footY]; dqe = [dq; dfootX; dfootY]; xe = [qe;dqe];
tau = [tau0; tau1];
kappa = [k 0 0 0; 0 0 0 0; 0 0 k 0; 0 0 0 0]; %spring stiffness
b = [beta 0 0 0; 0 0 0 0; 0 0 beta 0; 0 0 0 0]; %damping

%define position along stance leg; note that swing leg (leg 1) posn is a fn of
%stnace leg (leg 0) posn
%also find first time derivative of the position via chain rule, then get
%the jacobian of leg positions w/r.t. state vector q
alpha0 = int(q0,s,0,s) + tht0;
xsd0 = footX + d*D*cos(alpha0) - L*int(sin(alpha0),s,0,s);
ysd0 = footY + d*D*sin(alpha0) + L*int(cos(alpha0),s,0,s);
dxsd0 = diff(xsd0,q0)*dq0 + diff(xsd0,tht0)*dtht0;
dysd0 = diff(ysd0,q0)*dq0 + diff(ysd0,tht0)*dtht0;

%find hip positions for impact mapping
hipX = subs(xsd0,s,1);
hipY = subs(ysd0,s,1);

%find cartesian pos'n along swing leg; leg starts from hip and goes down
%hip pos'n is at the end of stance leg, so integrate along leg until you
%hit the end (s==1)
%with this formulation the swing leg will have positive values for tht0 &
%q0, the stance leg will have negative vals for tht1 and q1 b/c msrmnt from
%the hip
%because the orientation of the two legs are independent of one another
%both base angles tht need a group action applied to account for slope
alpha1 = int(q1,s,0,s) + tht1; %dont add alpha0 to this; the legs are connected @ hip by ideal joint
xsd1 = d*D*cos(alpha1) - L*int(sin(alpha1),s,0,s) + hipX; %leave this as just the foot position
ysd1 = d*D*sin(alpha1) + L*int(cos(alpha1),s,0,s) + hipY;
dxsd1 = diff(xsd1,q1)*dq1 + diff(xsd1,tht1)*dtht1;
dysd1 = diff(ysd1,q1)*dq1 + diff(ysd1,tht1)*dtht1;

%find swing foot vertical velocity for impact detection
dySw = jacobian(subs(ysd1,s,1),q);
dxSw = jacobian(subs(xsd1,s,1),q);
xSw = subs(xsd1,s,1);
ySw = subs(ysd1,s,1);

%get foot position

%find inertia matrix M
psd = jacobian([xsd0;ysd0;xsd1;ysd1],q);
fprintf("Nonextended Position Jacobian Found!\n");
M = m*(psd.'*psd);
M = int(M,d,-.5,.5);
M = int(M,s,0,1);
M = simplify(expand(M));
fprintf("Nonextended Inertia Matrix Found!\n")

%find extended inertia matrix Me for impact mapping
psde = jacobian([xsd0;ysd0;xsd1;ysd1;footX;footY],qe);
fprintf("Extended Position Jacobian Found!\n")
Me = m*(psde.'*psde);
Me = int(Me,d,-.5,.5);
Me = int(Me,s,0,1);
Me = simplify(expand(Me));
fprintf("Extended Inertia Matrix Found!\n");

%find E2 (swing foot pos'n deriv w/r.t. ext state vector) & dGamma (hip
%pos'n jacobian w/r.t. nonextended state vector) for impact mapping
E2 = jacobian(subs(ysd1,s,1),qe);
dGamma = [jacobian(hipX,q);jacobian(hipY,q)];

%find Gravity Vector & total energy expressions
%can be broken down into individual legs but would require reformulating M
%by leg to get KE/leg, same w PE
PE = int(int((m*g*(sin(phi)*xsd0 + cos(phi)*ysd0 + sin(phi)*xsd1 + cos(phi)*ysd1)),d,-0.5,0.5),s,0,1);
G = jacobian(PE,q).';
fprintf("Gravity Vector Found!\n")
KE = .5*(dq.'*M*dq);

E = KE + PE + .5*(q.'*kappa*q); 

%coriolis matrix - cristoffel symbols
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
C = simplify(C); %wrong btw - 10052022
fprintf("Coriolis Matrix Found!\n")

% find system lqr to check gains
% doesnt work with the extended formulation in symbolics
ddx = [dq;...
       M\([tau0;0;tau1;0] - C*dq - b*dq - kappa*q - G)];
A2 = jacobian(ddx,x);
B2 = jacobian(ddx,tau);
C2 = jacobian(x,x);
D2 = jacobian(x,tau);
fprintf("System Matrices Found!\n");

[D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars();
dq0 = 0; tht0 = 1e-15; dtht0 = 1e-15; tau0 = 0; q0 = 1e-30;
dq1 = 0; tht1 = 1e-15; dtht1 = 1e-15; tau1 = 0; q1 = 1e-30;
A2 = double(subs(A2));
B2 = double(subs(B2));
C2 = double(subs(C2));
D2 = double(subs(D2));
fprintf("System Matrix Substitution Done!\n");
sys = ss(A2,B2,C2,D2);
Q = eye(8);
fprintf("LQR Gains Found to BE:\n");
lqr(sys,Q,1)
fprintf("Calculations done, converting to file...\n")
% matlabFunction(M,C,G,E,PE,KE,"file","dynMats")%,"Optimize",false)
% matlabFunction(G,"file","groupActionFootForm")
% matlabFunction(Me,E2,dGamma,"file","impactMappingMatsFootForm")
% matlabFunction(dSwingFootX,dSwingFootY,swingFootX,swingFootY,"file","swingFootImpact")
% matlabFunction(E,"file","getEnergy")
