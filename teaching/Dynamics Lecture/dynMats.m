function [M,C,G] = dynMats(t,x)
%helper function to return the manipulator form matrices for use with ODE45
%simulates the physics of a double pendulum

%unpack inputs
tht1 = x(1); tht2 = x(2);
dtht1 = x(3); dtht2 = x(4);

t2 = cos(tht1);
t3 = sin(tht1);
t4 = tht1+tht2;
t5 = cos(t4);
t6 = sin(t4);
t7 = t5.*2.0;
t8 = t6.*2.0;
t11 = t6.*(9.81e+2./5.0e+1);
t9 = t2+t7;
t10 = t3+t8;
t12 = t7.*t9;
t13 = t8.*t10;
t14 = t6.*t9.*8.0;
t15 = t5.*t10.*8.0;
t16 = -t15;
t17 = t12+t13;
M = reshape([t2.^2+t3.^2+t9.^2+t10.^2,t17,t17,t5.^2.*4.0+t6.^2.*4.0],[2,2]);
if nargout > 1
    t18 = t14+t16;
    t19 = (dtht1.*t18)./4.0;
    C = reshape([dtht2.*t18.*(-1.0./4.0),t19,-t19+(dtht2.*(t5.*t10.*4.0-t6.*t9.*4.0))./2.0,0.0],[2,2]);
end
if nargout > 2
    G = [t3.*(9.81e+2./5.0e+1)+t11;t11];
end
