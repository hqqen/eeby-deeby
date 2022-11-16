function [M,C,G] = acrobotDyn(q)
%dynamics for acrobot problem
[m1, m2, l1, l2, I1, I2, g] = getSimParams();
q1 = q(1); q2 = q(2); dq1 = q(3); dq2 = q(4);

t2 = cos(q2);
t3 = sin(q1);
t4 = sin(q2);
t5 = q1+q2;
t6 = sin(t5);
t7 = (l1.*l2.*m2.*t2)./2.0;
t8 = I2+t7;
M = reshape([I1+I2+l1.^2.*m2+l1.*l2.*m2.*t2,t8,t8,I2],[2,2]);
if nargout > 1
    C = reshape([-dq2.*l1.*l2.*m2.*t4,(dq1.*l1.*l2.*m2.*t4)./2.0,dq2.*l1.*l2.*m2.*t4.*(-1.0./2.0),0.0],[2,2]);
end
if nargout > 2
    G = [g.*m2.*(l1.*t3+(l2.*t6)./2.0)+(g.*l1.*m1.*t3)./2.0,(g.*l2.*m2.*t6)./2.0];
end
