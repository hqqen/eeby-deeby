function [num] = costfun(x)
lambda = 0.125; 
k = (2.*pi)./lambda; 
f = [1;0.2];
xm = [15;195;285;105;15];
ym = [15;15;15;105;465];
xrec = [450.5;350.5];
yrec = [495.5;495.5];
N_a = length(xm);
N_s = length(xrec);
theta = atan(yrec./xrec);
a = x(1:N_a,1);
alpha = x(N_a+1:2*N_a,1);

d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-[xrec(i,:);yrec(i,:)]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end

u = zeros(N_a,N_s);
v = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        u(m,i) = 1/d(m,i)*cos(alpha(m,:)+zeta(m,i));
        v(m,i) = 1/d(m,i)*sin(alpha(m,:)+zeta(m,i));
    end
end
for i=1:N_s
    den = 0;
    for m=1:N_a
        den = den + a(m,:)*(u(m,i)+1i*v(m,i));
    end
    den = abs(den);
    num(i,:) = den - f(i,:);
end
end

