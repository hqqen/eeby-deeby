clc
clear all

%%
lambda = 7.5; 
k = (2.*pi)./lambda; 
xm = [0;2.65;5.3;7.95;10.6;13.25]+30;
ym = [0;2.6;5.3;7.9;10.6;13.2]+220;
xrec = 1:280;
yrec = 2:250;
lx = length(xrec);
ly = length(yrec);
for i=1:lx
    xv(i,1) = xrec(1,i);
    yv(i,1) = 1;
end
for i=1:ly
    xv(lx+i,1) = 280;
    yv(lx+i,1) = yrec(1,i);
end
xrec = xv;
yrec = yv;
N_a = length(xm);
N_s = length(xrec);
% scatter(xrec,yrec)
rho = sqrt(xrec.^2+yrec.^2);
theta = atan(yrec./xrec);
load('x_sensor.mat')
load('y_sensor.mat')
for i=1:N_s
    pos_x(i,1) = find(xv == xrec(i,1));
    pos_y(i,1) = find(yv == yrec(i,1));
end
for m=1:N_a
    if m==1
        load('amat_sensor_0.mat');
    elseif m==2
        load('amat_sensor_1.mat');
    elseif m==3
        load('amat_sensor_2.mat');
    elseif m==4
        load('amat_sensor_3.mat');
    elseif m==5
        load('amat_sensor_4.mat');
    elseif m==6
        load('amat_sensor_5.mat');
    end
    for i=1:N_s
        gamma_ch(m,i) = 1/0.35*amat_fdtd(pos_x(i,1),pos_y(i,1));
    end
end
mu = 2;
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-[xrec(i,:);yrec(i,:)]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end
load('a.mat');
load('alpha.mat');
a = as;
alpha = alphas;
for m=1:N_a
    for i=1:N_s
        u_ch(m,i) = gamma_ch(m,i)*cos(alpha(m,1)+zeta(m,i));
        v_ch(m,i) = gamma_ch(m,i)*sin(alpha(m,1)+zeta(m,i));
    end
end
for i=1:N_s
    den_ch = 0;
    for m=1:N_a
        den_ch = den_ch + a(m,1)*(u_ch(m,i)+1i*v_ch(m,i));
    end
    den_ch = abs(den_ch);
    af(i,:) = den_ch;
end
af_db = 20*log10(af/max(af));
plot3(xrec,yrec,af_db, 'LineWidth', 3)
grid on
hold on
scatter3(xm,ym,-40*ones(6,1), 'LineWidth', 3)
scatter3([220;240],[30;10],[-40;-40],'d', 'LineWidth', 3)
scatter3([270],[30],[-40],'+', 'LineWidth', 3)
xlabel('x')
ylabel('y')
set(gca, 'LineWidth', 5, 'FontSize', 35)




