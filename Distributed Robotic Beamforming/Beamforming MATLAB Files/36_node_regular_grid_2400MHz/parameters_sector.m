clc
% clear all

%%
lambda = 0.125; 
k = (2.*pi)./lambda; 
xm = [105;195;196;285;285];
ym = [195;200;285;195;285];
s = 3;
f = [1;0.1*ones(s,1)];
w = 1./f;
w = w/sum(w);
r = 150;
theta = pi./12;
theta_sec = linspace(pi./3,1/2.*pi,s);
theta = [theta;theta_sec'];
xrec = 225 + r*cos(theta);
yrec = 225 + r*sin(theta);
xrec = round(xrec) + 0.5;
yrec = round(yrec) + 0.5;
N_a = length(xm);
N_s = length(xrec);
rho = sqrt(xrec.^2+yrec.^2);
theta = atan(yrec./xrec);
xv = (0.5:1:499.5)';
yv = xv;
for i=1:N_s
    pos_x(i,1) = find(xv == xrec(i,1));
    pos_y(i,1) = find(yv == yrec(i,1));
end
for m=1:N_a
    if m==1
        load('amat_105_195.mat');load('pmat_105_195.mat');
    end
    if m==2
        load('amat_195_200.mat');load('pmat_195_200.mat');
    end
    if m==3
        load('amat_196_285.mat');load('pmat_196_285.mat');
    end
    if m==4
        load('amat_285_195.mat');load('pmat_285_195.mat');
    end
    if m==5
        load('amat_285_285.mat');load('pmat_285_285.mat');
    end
    amat = amat/max(amat);
    amat = reshape(amat,500,500);
    pmat = reshape(pmat,500,500);
    for i=1:N_s
        gamma_amp(m,i) = 1/1*amat(pos_x(i,1),pos_y(i,1));
        gamma_ph(m,i) = 1/1*pmat(pos_x(i,1),pos_y(i,1));
    end
end
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
scatter(xm,ym,'LineWidth',5)
hold on
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
scatter(rho(1).*cos(theta(1)),rho(1).*sin(theta(1)),'LineWidth',5)
scatter(rho(2:end).*cos(theta(2:end)),rho(2:end).*sin(theta(2:end)),'LineWidth',5)
xlabel('position x')
ylabel('position y')
axis('square')
xlim([50 400])
ylim([50 400])
% legend('T_x','R_x','Adversary')