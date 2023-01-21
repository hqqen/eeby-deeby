syms a alpha x y [1 4] real
syms gamma d k tht real

AF = sum( ( (a.'*gamma.'/d).'*(exp(1i.'*(alpha + k.'*x.'*cos(tht) + k.'*y.'*sin(tht) + k.'*d))) ) );
ddAF = hessian(AF, [a alpha])