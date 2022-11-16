function [M,C,G,T,V] = getDynamics(q,dq)
%helper function to return MPFL normal form compliant dynamics along with
%kinetic and potential energy
%takes state q and first time deriv dq
%returns state matrices
tht1 = q(1); tht2 = q(2); dtht1 = dq(1); dtht2 = dq(2);
[I1,I2,m1,m2,a1,a2,k1,k2] = getVars();

M = reshape([I1+I2,I2,I2,I2],[2,2]);
if nargout > 1
    C = reshape([0.0,0.0,0.0,0.0],[2,2]);
end
if nargout > 2
    t2 = cos(tht1);
    t3 = sin(tht1);
    t4 = tht1+tht2;
    t5 = sin(t4);
    G = [k1.*tht1+m2.*(a1.*t3+a2.*t5)+a1.*m1.*t3;k2.*tht2+a2.*m2.*t5];
end
if nargout > 3
    T = (I2.*(dtht1+dtht2).^2)./2.0+(I1.*dtht1.^2)./2.0;
end
if nargout > 4
    V = (k1.*tht1.^2)./2.0+(k2.*tht2.^2)./2.0-m2.*(a1.*t2+a2.*cos(t4))-a1.*m1.*t2;
end
