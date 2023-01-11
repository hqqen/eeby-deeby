function [M,C,G,P] = dynMats(q,dq)
%helper function to numerically calculate dynamics matrices M (inertia), C
%(coriolis) and G (gravity) given state vector q and its derivatives

q0 = q(1); tht = q(2);
dq0 = dq(1); dtht = dq(2);

% fprintf("q:\n"); disp(q);
% fprintf("\ndq:\n"); disp(dq);
% fprintf("\nqo:\n"); disp(q0);
% fprintf("\ntht:\n"); disp(tht);
% fprintf("\ndq0:\n"); disp(dq0);
% fprintf("\ndtht:\n"); disp(dtht);

[D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars();

t2 = cos(q0);
t3 = sin(q0);
t4 = phi+tht;
t5 = D.^2;
t6 = L.^2;
t7 = q0.^2;
t8 = q0.^3;
t9 = q0.^5;
t19 = 1.0./q0.^7;
t10 = cos(t4);
t11 = q0.*t3;
t12 = q0+t4;
t13 = sin(t4);
t14 = 1.0./t7;
t15 = 1.0./t8;
t17 = 1.0./t9;
t24 = t6.*7.2e+1;
t25 = t5.*t9;
t26 = m.*q0.*t6.*2.0;
t27 = m.*t3.*t6.*2.0;
t30 = (m.*t5)./2.4e+1;
t31 = t6.*t8.*1.2e+1;
t34 = t3.*t6.*1.44e+2;
t37 = t2.*t6.*-7.2e+1;
t16 = t14.^2;
t18 = t14.^3;
t20 = cos(t12);
t21 = sin(t12);
t22 = q0.*t13;
t23 = -t10;
t28 = q0.*t24;
t29 = -t26;
t32 = t2.*t24;
t33 = t2+t11-1.0;
t36 = t11.*t24;
t38 = -t34;
t39 = m.*t2.*t6.*t15;
t40 = m.*t6.*t17.*1.2e+1;
t41 = (m.*t6.*t14)./2.0;
t43 = m.*t2.*t6.*t17.*4.8e+1;
t35 = t2.*t28;
t42 = -t39;
t44 = m.*t3.*t6.*t16.*1.0e+1;
t45 = t27+t29;
t46 = m.*t6.*t16.*t33;
t47 = t20+t22+t23;
t49 = m.*t6.*t19.*t33.*1.2e+2;
t48 = -t46;
t50 = -t49;
t52 = t25+t28+t31+t35+t38;
t51 = t30+t41+t48;
M = reshape([(m.*t17.*t52)./3.6e+1,t51,t51,(m.*t5)./1.2e+1+t15.*(t26-t27)],[2,2]);
if nargout > 1
    t53 = t40+t42+t43+t44+t50;
    C = reshape([-dq0.*((dtht.*((dq0.*t53)./2.0-(dtht.*(t39.*2.0+t17.*(m.*t6.*2.0-m.*t2.*t6.*2.0).*3.6e+1-t18.*(t26-t27).*6.0e+1-m.*t3.*t6.*t16.*1.8e+1))./2.0))./2.0+(dq0.*((dtht.*t53)./2.0+dq0.*m.*t18.*(t28+q0.*t37+t5.*t8.*2.0e+1).*(5.0./2.4e+1)-dq0.*m.*t19.*(t24+t37+t6.*t7.*3.6e+1-t6.*t11.*7.2e+1+(t5.*5.0)./t16).*(5.0./4.0)+dq0.*m.*t16.^2.*t52.*(3.5e+1./1.2e+1)-(dq0.*m.*t17.*(t24+t36+t37+t5.*t7.*6.0e+1))./7.2e+1))./2.0),0.0,0.0,0.0],[2,2]);
end
if nargout > 2
    G = [L.*g.*m.*t15.*t47.*2.0-L.*g.*m.*t14.*(t13-t21);-L.*g.*m.*t14.*(t13-t21+q0.*t10)];
end
if nargout > 3
    P = -L.*g.*m.*t14.*t47;
end
