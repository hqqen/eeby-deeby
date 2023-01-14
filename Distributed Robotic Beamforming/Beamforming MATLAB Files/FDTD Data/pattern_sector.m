clc
clear all

%%
lambda = 7.5; 
k = (2.*pi)./lambda; 
xm = [0;2.65;5.3;7.95;10.6;13.25;18.55]+30;
ym = [0;2.6;5.3;7.9;10.6;13.2;18.5]+220;
load('x_sensor.mat')
load('y_sensor.mat')
r = (max(yv)-min(yv));
s = 100;
theta = linspace(-pi/2,0,s);
xsec = min(xv)+r*cos(theta);
ysec = max(yv)+r*sin(theta);
xsec = round(xsec);
ysec = round(ysec);
xsec(1) = 1;
ysec(1:3) = ones(3,1);
xsec(end) = 264;
ysec(end) = max(yv);
xrec = xsec';
yrec = ysec';
N_a = length(xm);
N_s = length(xrec);
% scatter(xrec,yrec)
rho = sqrt(xrec.^2+yrec.^2);
theta = atan(yrec./xrec);
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
        gamma_ch(m,i) = 1*amat_fdtd(pos_x(i,1),pos_y(i,1));
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
load('a_sector.mat')
load('alpha_sector.mat')
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
af = [af;1.6;3.3;1.8];
af_db = 20*log10(af/max(af));

%%
plot3(xrec,yrec,af_db(1:s), 'LineWidth', 3)
grid on
hold on
scatter3(xm,ym,-40*ones(7,1), 'LineWidth', 3)
scatter3([220;240],[30;10],[af_db(s+1);af_db(s+2)],'d', 'LineWidth', 3)
scatter3([270],[30],[af_db(s+3)],'+', 'LineWidth', 3)
xlabel('x')
ylabel('y')
set(gca, 'LineWidth', 5, 'FontSize', 35)