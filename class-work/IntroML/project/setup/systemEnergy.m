function [H,L,T,U] = systemEnergy(q)
%system energy for acrobot problem
[m1, m2, l1, l2, I1, I2, g] = getSimParams();
q1 = q(1); q2 = q(2); dq1 = q(3); dq2 = q(4);

t2 = cos(q1);
t3 = cos(q2);
t4 = q1+q2;
t5 = dq1.^2;
t6 = dq2.^2;
t7 = l1.^2;
t8 = l1.*t2;
t9 = cos(t4);
t10 = (I1.*t5)./2.0;
t11 = (I2.*t5)./2.0;
t12 = (I2.*t6)./2.0;
t15 = (l1.*l2.*m2.*t3)./2.0;
t18 = (m2.*t5.*t7)./2.0;
t13 = (l2.*t9)./2.0;
t14 = (g.*m1.*t8)./2.0;
t17 = I2+t15;
t19 = t5.*t15;
t16 = -t14;
t20 = dq1.*dq2.*t17;
t21 = t8+t13;
t22 = g.*m2.*t21;
t23 = -t22;
H = t10+t11+t12+t16+t18+t19+t20+t23;
if nargout > 1
    L = t10+t11+t12+t14+t18+t19+t20+t22;
end
if nargout > 2
    T = t10+t11+t12+t18+t19+t20;
end
if nargout > 3
    U = t16+t23;
end
