function G = getGroupAction(q0,q1,tht0,tht1)
%getGroupAction
%    G = getGroupAction(Q0,Q1,THT0,THT1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Jan-2023 13:03:51

t2 = cos(tht0);
t3 = cos(tht1);
t4 = sin(tht0);
t5 = sin(tht1);
t6 = q0+tht0;
t7 = q1+tht1;
t8 = q1.^2;
t9 = q1.^3;
t13 = 1.0./q0.^2;
t10 = cos(t6);
t11 = cos(t7);
t12 = sin(t6);
t14 = 1.0./t9;
t15 = 1.0./t8.^2;
t16 = t3.*(9.8e+1./5.0);
t18 = q1.*t5.*(9.8e+1./5.0);
t20 = t2.*t9.*(4.9e+1./5.0);
t22 = t4.*t9.*(9.8e+1./5.0);
t17 = -t16;
t19 = t11.*(9.8e+1./5.0);
t21 = -t20;
t23 = t9.*t10.*(4.9e+1./5.0);
t24 = t9.*t12.*(4.9e+1./5.0);
t25 = -t24;
t26 = t22+t25;
t27 = q0.*t26;
t28 = t21+t23+t27;
G = [t13.*t14.*(-t22+q0.*t23+t9.*t12.*(9.8e+1./5.0))+1.0./q0.^3.*t14.*t28.*2.0,t13.*t14.*(t24-t4.*t9.*(4.9e+1./5.0)+q0.*(t23-t2.*t9.*(9.8e+1./5.0))),t14.*(t17+t18+t19)+t15.*(t5.*(9.8e+1./5.0)-sin(t7).*(9.8e+1./5.0)+q1.*t16-t5.*t8.*(4.9e+1./5.0)).*3.0+t13.*t15.*t28.*3.0-t13.*t14.*(t2.*t8.*(-1.47e+2./5.0)+t8.*t10.*(1.47e+2./5.0)+q0.*(t4.*t8.*(2.94e+2./5.0)-t8.*t12.*(1.47e+2./5.0))),t14.*(t17+t18+t19+t3.*t8.*(4.9e+1./5.0))];
end
