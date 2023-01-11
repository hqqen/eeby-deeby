clc
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,T_out);ym(m,T_out)]-rho*[cos(theta(i,:));sin(theta(i,:))]);
        zeta(m,i) = k*(xm(m,T_out)*cos(theta(i,:)) + ym(m,T_out)*sin(theta(i,:)) + d(m,i));
    end
end
den = zeros(N_s,1);
u = zeros(N_a,N_s);
v = zeros(N_a,N_s);
t = T-1;
x = [a(:,t);alpha(:,t)];
% a = x(1:N_a,:);
% alpha = x(N_a+1:2*N_a,:);
for m=1:N_a
    for i=1:N_s
        u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i));
        v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
    end
end
for i=1:N_s
    for m=1:N_a
        den(i,:) = den(i,:) + x(m,1)*(u(m,i)+1i*v(m,i));
    end
end
den = abs(den);
error_f = norm(den-f)/N_s;
den_db = 10*log10(den);

% h = polar(theta,den,'r');
% h.LineWidth = 2;
% hold on

plot(den_db,'b','LineWidth',5)
hold on
grid on

% xlabel('iteration t')
% ylabel('||g(t)||')
% set(gca, 'LineWidth', 5, 'FontSize', 35)