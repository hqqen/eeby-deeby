clc
clear all

%%
lambda = 7.5; 
k = (2.*pi)./lambda; 
f = [100;100];
w = 1./f;
w = w/sum(w);
xm = [0;2.65;5.3;7.95;10.6;13.25;18.55]+30;
ym = [0;2.6;5.3;7.9;10.6;13.2;18.5]+220;
% xm = [0;7.95;13.25;18.55];
% ym = [0;7.9;13.2;18.5];
xrec = [220;240;260];
yrec = [10;30;10];
N_a = length(xm);
N_s = length(xrec);
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
    elseif m==7
        load('amat_sensor_7.mat');
    end
    for i=1:N_s
        gamma_ch(m,i) = 1/0.4*amat_fdtd(pos_x(i,1),pos_y(i,1));
    end
end
% gamma_ch = ones(N_a,N_s);
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
e1 = gamma_ch.*d
