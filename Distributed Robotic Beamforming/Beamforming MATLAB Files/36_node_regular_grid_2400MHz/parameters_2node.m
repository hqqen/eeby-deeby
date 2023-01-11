clc
% clear all

%%
lambda = 7.5; 
k = (2.*pi)./lambda; 
f = 50;
w = 1;
xm = [0;8];
ym = [35;35];
xrec = [0];
yrec = [0];
N_a = length(xm);
N_s = length(xrec);
rho = sqrt(xrec.^2+yrec.^2);
theta = pi;
mu = 2;

%%
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-[xrec(i,:);yrec(i,:)]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end

%%
scatter(xm,ym,100,'filled')
hold on
scatter(xrec(1),yrec(1),100,'filled','d','g')
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
xlabel('x-coordinate')
ylabel('y-coordinate')
legend('T_x','R_x')
        
