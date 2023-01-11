clear af_sec
s_sec = 50;
rho_sec = [50*lambda/2*ones(s_sec,1)];
theta_sec = linspace(0,2.*pi,s_sec);
theta_sec = theta_sec';
gamma_ch = normrnd(1,0.1,N_a,s_sec);
as = a(:,q);
alphas = alpha(:,q);
d = zeros(N_a,s_sec);
zeta = zeros(N_a,s_sec);
for m=1:N_a
    for i=1:s_sec
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho_sec(i,:).*[cos(theta_sec(i,:));sin(theta_sec(i,:))]);
        zeta(m,i) = k*(xm(m,:)*cos(theta_sec(i,:)) + ym(m,:)*sin(theta_sec(i,:)) + d(m,i));
    end
end
for m=1:N_a
        for i=1:s_sec
            u_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*cos(alphas(m,1)+zeta(m,i));
            v_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*sin(alphas(m,1)+zeta(m,i));
        end
    end
for i=1:s_sec
    den_ch = 0;
    for m=1:N_a
        den_ch = den_ch + as(m,1)*(u_ch(m,i)+1i*v_ch(m,i));
    end
    den_ch = abs(den_ch);
    af_sec(i,:) = den_ch;
end
x_sec = rho_sec.*cos(theta_sec);
y_sec = rho_sec.*sin(theta_sec);
close all
plot3(x_sec,y_sec,20*log10(af_sec/max([af_sec;af(1:2,q)])), 'LineWidth', 5)
hold on
xrec = rho(1:2).*cos(theta(1:2));
yrec = rho(1:2).*sin(theta(1:2));
scatter3(xrec,yrec,20*log10(af(1:2,q)/max([af_sec;af(1:2,q)])),'d', 'LineWidth', 5)
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
xlabel('position x')
ylabel('position y')
zlabel('R_x power(dB)')