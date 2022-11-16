function [p1,p2] = getJointPosition(q)
%joint positions for acrobot problem
%mostly deprecated
[~, ~, l1, l2, ~, ~, ~] = getSimParams();
q1 = q(1); q2 = q(2); dq1 = q(3); dq2 = q(4);

t2 = cos(q1);
t3 = sin(q1);
t4 = q1+q2;
t5 = l1.*t2;
t6 = l1.*t3;
t7 = -t5;
p1 = [t6;t7];
if nargout > 1
    p2 = [t6+l2.*sin(t4);t7-l2.*cos(t4)];
end
