clc
clear all

%%
lambda = 7.5; 
k = (2.*pi)./lambda; 
for i=1:80
    if rand<0.5
        xm(i,:) = (mod(i,7)-3)*1;
        ym(i,:) = (ceil(i/7)-4)*1;
    end
end
ind = find(xm);
xm = xm(ind);
ym = ym(ind);
xm = xm*lambda/2;
ym = ym*lambda/2;
s = 20;
f = [5;1;2*ones(s,1)]/10;
rho = [50*lambda/2*ones(2,1);15*lambda/2*ones(s,1)];
theta = [pi./6;pi./4];
theta_sec = linspace(0,2.*pi,s);
theta = [theta;theta_sec'];
N_a = length(xm);
N_s = length(rho);
mu = 2;

%%
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho(i,:).*[cos(theta(i,:));sin(theta(i,:))]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end

%%
scatter(xm,ym,'LineWidth',5)
hold on
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
scatter(rho(1:2).*cos(theta(1:2)),rho(1:2).*sin(theta(1:2)),'LineWidth',5)
scatter(rho(3:end).*cos(theta(3:end)),rho(3:end).*sin(theta(3:end)),'LineWidth',5)
xlabel('position x')
ylabel('position y')
axis('square')
ylim([-100 200])