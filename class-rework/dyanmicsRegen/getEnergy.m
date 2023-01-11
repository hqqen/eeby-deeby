function [T,P] = getEnergy(dq0,dq1,dtht0,dtht1,q0,q1,tht0,tht1)
%getEnergy
%    [T,P] = getEnergy(DQ0,DQ1,DTHT0,DTHT1,Q0,Q1,THT0,THT1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Jan-2023 11:59:49

t2 = cos(q0);
t3 = cos(q1);
t4 = sin(q0);
t5 = sin(q1);
t6 = q0+tht0;
t7 = q1.*2.0;
t8 = q0.^2;
t9 = q1.^2;
t12 = -q1;
t13 = 1.0./q0;
t17 = 1.0./q1.^3;
t19 = -tht0;
t20 = -tht1;
t10 = t9.^2;
t11 = q0.*t4;
t14 = 1.0./t8;
t15 = 1.0./t9;
t21 = t20+tht0;
t23 = t6+t20;
t24 = q1+t19+tht1;
t26 = t3.*2.4e+3;
t27 = t2-3.0./2.0;
t28 = t9.*1.2e+3;
t34 = q1.*t5.*2.4e+3;
t16 = t14.^2;
t18 = 1.0./t10;
t22 = cos(t21);
t25 = sin(t21);
t29 = -t26;
t30 = cos(t23);
t31 = cos(t24);
t32 = sin(t23);
t33 = sin(t24);
t37 = -t34;
t49 = t8.*t27;
t53 = t12+t23;
t35 = t22.*2.0;
t36 = t25.*2.0;
t39 = q1.*t25;
t40 = q1.*t30;
t42 = q1.*t32;
t43 = q1.*t33;
t45 = t30.*2.0;
t46 = t31.*2.0;
t47 = t32.*2.0;
t48 = t33.*2.0;
t50 = -t30;
t51 = -t31;
t52 = t12.*t22;
t54 = -t32;
t56 = t12.*t25;
t57 = cos(t53);
t58 = t12.*t31;
t59 = sin(t53);
t60 = t12.*t32;
t71 = t2+t11+t49-1.0;
t75 = t10+t28+t29+t37+2.4e+3;
t44 = -t35;
t55 = -t47;
t61 = q1.*t57;
t62 = q1.*t59;
t63 = t57.*2.0;
t64 = t59.*2.0;
t66 = t40./1.2e+3;
t67 = t42./1.2e+3;
t73 = (t9.*t59)./1.2e+3;
t74 = t16.*t71;
t77 = t40+t54+t59;
t65 = -t63;
t68 = -t67;
t69 = t61./1.2e+3;
t70 = t62./1.2e+3;
t76 = t74-1.0./8.0e+2;
t78 = q0.*t77;
t84 = t25+t33+t52+t77;
t72 = -t69;
t79 = t68+t70;
t81 = t42+t45+t62+t65;
t86 = t13.*t15.*t84;
t88 = t22+t39+t50+t51+t57+t60+t78;
t80 = t15.*t79;
t82 = q0.*t81;
t83 = t66+t72+t73;
t87 = t43+t44+t46+t56+t81;
t90 = t14.*t15.*t88;
t85 = t17.*t83;
t89 = t13.*t17.*t87;
t91 = t80+t86;
t92 = t36+t40+t48+t52+t55+t58+t61+t64+t82;
t93 = t80+t90;
t94 = t14.*t17.*t92;
t95 = t85+t89;
t96 = t85+t94;
et1 = -dq1.*((dq0.*t96)./2.0+(dtht0.*t95)./2.0-(dtht1.*t18.*t75)./4.8e+3-(dq1.*1.0./q1.^5.*(t5.*-4.0+t7+t3.*t7+q1.^5./3.6e+3+1.0./(t17.*3.0)))./2.0)-dtht1.*((dq0.*t93)./2.0+(dtht0.*t91)./2.0+(dtht1.*(t17.*(t5.*2.0-t7)-8.333333333333333e-4))./2.0-(dq1.*t18.*t75)./4.8e+3)-dtht0.*((dq0.*t76)./2.0+(dq1.*t95)./2.0+(dtht1.*t91)./2.0+(dtht0.*(t13.^3.*(t4.*2.0+q0.*(t2.*2.0-4.0))-1.0./6.0e+2))./2.0);
et2 = -dq0.*((dq1.*t96)./2.0+(dtht0.*t76)./2.0+(dtht1.*t93)./2.0-(dq0.*(t13.^5.*(q0.*4.0-t4.*4.0-t4.*t8.*2.0+1.0./t13.^3.*(4.0./3.0))+1.0./9.0e+2))./2.0);
T = et1+et2;
if nargout > 1
    P = -t15.*(cos(q1+tht1).*(4.9e+1./5.0)-cos(tht1).*(4.9e+1./5.0)+q1.*sin(tht1).*(4.9e+1./5.0))+t14.*t15.*(q0.*(t9.*sin(t6).*(4.9e+1./5.0)-t9.*sin(tht0).*(9.8e+1./5.0))-t9.*cos(t6).*(4.9e+1./5.0)+t9.*cos(tht0).*(4.9e+1./5.0));
end
end
