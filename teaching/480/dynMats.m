function [M,C,G] = dynMats(dtht1,dtht2,tht1,tht2)
%dynMats
%    [M,C,G] = dynMats(DTHT1,DTHT2,THT1,THT2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    06-Dec-2022 19:10:58

t2 = cos(tht1);
t3 = sin(tht1);
t4 = sin(tht2);
t5 = tht1+tht2;
t6 = dtht1.*2.0;
t7 = cos(t5);
t8 = sin(t5);
t9 = dtht2+t6;
t10 = t7.*2.0;
t11 = t8.*2.0;
t14 = t8.*(9.81e+2./5.0e+1);
t12 = t2+t10;
t13 = t3+t11;
t15 = t10.*t12;
t16 = t11.*t13;
t17 = t15+t16;
M = reshape([t2.^2+t3.^2+t12.^2+t13.^2,t17,t17,t7.^2.*4.0+t8.^2.*4.0],[2,2]);
if nargout > 1
    C = reshape([0.0,t4.*t9,t4.*t9.*-2.0,-dtht1.*t4],[2,2]);
end
if nargout > 2
    G = [t3.*(9.81e+2./5.0e+1)+t14;t14];
end
