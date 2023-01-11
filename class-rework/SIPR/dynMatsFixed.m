function [M,C,G,E,PE,KE] = dynMatsFixed(q,dq)
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
t7 = dq0.^2;
t8 = dtht.^2;
t9 = q0.*2.0;
t10 = q0.^2;
t11 = q0.^3;
t21 = 1.0./q0.^5;
t25 = q0./2.0;
t12 = t3.*3.0;
t13 = q0.*t2;
t14 = cos(t4);
t15 = q0.*t3;
t16 = q0+t4;
t17 = sin(t4);
t18 = 1.0./t10;
t19 = 1.0./t11;
t26 = sin(t25);
t30 = (m.*t5)./2.4e+1;
t38 = (m.*t5.*t7)./7.2e+1;
t43 = m.*t3.*t6.*t7.*t21.*2.0;
t20 = t18.^2;
t22 = cos(t16);
t23 = sin(t16);
t24 = -t12;
t27 = q0.*t17;
t28 = -t14;
t29 = t26.^2;
t31 = t2+t15-1.0;
t32 = dq0.*dtht.*t30;
t33 = m.*t6.*t8.*t18;
t34 = t8.*t30;
t35 = (m.*t6.*t18)./2.0;
t39 = dq0.*dtht.*m.*t3.*t6.*t19;
t41 = m.*t3.*t6.*t8.*t19;
t42 = (m.*t6.*t7.*t18)./6.0;
t45 = -t43;
t36 = m.*t6.*t7.*t20.*2.0;
t37 = t9+t13+t24;
t40 = dq0.*dtht.*t35;
t44 = -t39;
t46 = -t41;
t47 = m.*t6.*t20.*t31;
t48 = t22+t27+t28;
t50 = dq0.*dtht.*m.*t6.*t20.*t29.*2.0;
t52 = m.*t6.*t7.*t20.*t29.*-2.0;
t49 = -t47;
t51 = t29.*t36;
t53 = L.*g.*m.*t18.*t48;
t54 = -t53;
t55 = t30+t35+t49;
M = reshape([(m.*t5)./3.6e+1+t21.*(m.*t3.*t6.*-4.0+(m.*t6.*t11)./3.0+(m.*q0.*(t6.*7.2e+1+t2.*t6.*7.2e+1))./3.6e+1),t55,t55,(m.*t5)./1.2e+1-t19.*(m.*t3.*t6.*2.0-m.*t6.*t9)],[2,2]);
if nargout > 1
    C = reshape([dq0.*m.*t6.*t18.^3.*(q0.*1.2e+1-t3.*3.0e+1+t11+t13.*1.8e+1+t10.*t12).*(-1.0./3.0),-m.*t6.*t21.*(dq0.*4.0-dq0.*t2.*4.0+dq0.*t10-dq0.*t15.*4.0-dtht.*t15.*3.0+dtht.*q0.*t9+dq0.*t2.*t10+dtht.*t2.*t10),dtht.*m.*t6.*t20.*t37,-dq0.*m.*t6.*t20.*t37],[2,2]);
end
if nargout > 2
    G = [L.*g.*m.*t19.*t48.*2.0-L.*g.*m.*t18.*(t17-t23);-L.*g.*m.*t18.*(t17-t23+q0.*t14)];
end
if nargout > 3
    E = t32+t33+t34+t36+t38+t40+t42+t44+t45+t46+t50+t52+t54+(k.*t10)./2.0;
end
if nargout > 4
    PE = t54;
end
if nargout > 5
    KE = t32+t33+t34+t36+t38+t40+t42+t44+t45+t46+t50+t52;
end
