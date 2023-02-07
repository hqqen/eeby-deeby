clear all
syms a gamma mu tht x y rho k alpha real
d = sqrt((x - rho*cos(tht))^2 + (y - rho*sin(tht))^2);
zeta = k*(x*cos(tht) + y*sin(tht) + d);
du = a*gamma*[-(mu*cos(alpha + zeta))*(2*x - 2*rho*cos(tht))/(4*d^(2+mu/2)) - sin(alpha + zeta)*(k*cos(tht) + k*(2*x - 2*rho*cos(tht))/(2*d))*d^(-mu/2); ...
              -(mu*cos(alpha + zeta))*(2*y - 2*rho*sin(tht))/(4*d^(2+mu/2)) - sin(alpha + zeta)*(k*sin(tht) + k*(2*x - 2*rho*sin(tht))/(2*d))*d^(-mu/2)];

ddu = jacobian(du,[x,y]);
dduOut = latex(ddu);

dv = a*gamma*[-(mu*sin(alpha + zeta))*(2*x - 2*rho*cos(tht))/(4*d^(2+mu/2)) + cos(alpha + zeta)*(k*cos(tht) + k*(2*x - 2*rho*cos(tht))/(2*d))*d^(-mu/2); ...
              -(mu*sin(alpha + zeta))*(2*y - 2*rho*sin(tht))/(4*d^(2+mu/2)) + cos(alpha + zeta)*(k*sin(tht) + k*(2*x - 2*rho*sin(tht))/(2*d))*d^(-mu/2)];

ddv = jacobian(dv,[x,y]);
ddvOut = latex(ddv);