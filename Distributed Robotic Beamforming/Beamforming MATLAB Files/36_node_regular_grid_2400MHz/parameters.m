clc
% clear all

%%
lambda = 0.125; 
k = (2.*pi)./lambda; 
f = [1;0.25];
w = 1./f;
w = w/sum(w);
% xm = [15;195;285;105;15;15];
% ym = [15;15;15;105;195;465];
xm = [15;195;285;105;15];
ym = [15;15;15;105;465];
% xrec = [450.5;350.5];
% yrec = [495.5;495.5];
% xrec = [375.5;425.5];
% yrec = [495.5;425.5];
% xrec = [425.5;375.5];
% yrec = [425.5;495.5];
xrec = [375.5;375.5];
yrec = [375.5;495.5];
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
        load('amat_15_15.mat');load('pmat_15_15.mat');
    end
    if m==2
        load('amat_195_15.mat');load('pmat_195_15.mat');
    end
    if m==3
        load('amat_285_15.mat');load('pmat_285_15.mat');
    end
    if m==4
        load('amat_105_105.mat');load('pmat_105_105.mat');
    end
%     if m==5
%         load('amat_15_195.mat');load('pmat_15_195.mat');
%     end
    if m==5
        load('amat_15_465.mat');load('pmat_15_465.mat');
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
scatter(xm,ym,100,'filled')
hold on
scatter(xrec(1),yrec(1),100,'filled','d','g')
scatter(xrec(2),yrec(2),100,'filled','d','r')
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
xlabel('x-coordinate')
ylabel('y-coordinate')
legend('T_x','R_x','Null R_x')
xlim([0 500])
ylim([0 500])
        
