clc
clear all
close all

%%
lambda = 7.5;                       % inter agent spacing?
k = (2.*pi)./lambda;                % wave number
for i=1:60                          % randomize inital positions of agents
    if rand<0.5
        xm(i,:) = (mod(i,7)-3)*1;
        ym(i,:) = (ceil(i/7)-4)*1;
    end
end
ind = find(xm);                     % filter out any agents initialized at the origin
xm = xm(ind);
ym = ym(ind);
xm = xm*lambda/2;
ym = ym*lambda/2;
s = 20;                             % numer of sectors to look at
% f = [5;1*ones(s,1)]/10;             % set magnitude of des'd beampower at each receiver
rho = [50*lambda/2*ones(1,1);50*lambda/2*ones(s/2,1);60*lambda/2*ones(s/2,1)]; % get distance from transmitter to receiver
theta = pi./15;                         % angle of allied reciever
theta_sec = linspace(pi./3,3*pi./3,s);  % divide workspace into 30deg sectors
theta = [theta;theta_sec'];             % make list of all receiver posns
% theta = 0;                         % angle of allied reciever
% theta_sec = linspace(pi./3,2.5*pi./3,s);  % divide workspace into 30deg sectors
% theta = [theta;theta_sec'];             % make list of all receiver posns

% % Parameters of desired beam
% f = 40e6;                                     % Frequency in Hz                                                                        
% N_d = 5;                                       % Number of elements of array  
% lambda = 3e8/f;                               % wave length
% d = lambda/4;                                 % Inter-element spacing
% k = (2*pi)./lambda;                           % wave number
% rm_d = [xm,ym];
% % channel parameters
% gamma = 1;                                   
% mu = 2;
% nu = 1;  
% rho = 1.5*d; 
% % Array factor (AFd) of desired beam
% AFd = GetArrayFactor(N_d, rm_d, Am_d, alpha_d, gamma, mu, nu, rho, theta); % call provided code to find AF
% figure();
% h = polar(theta,abs(AFd),'r'); % plot
% set(h,'LineWidth',2);
% box on
% hold on


f = squareFourier(5,pi,.1,5,theta);             % set magnitude of des'd beampower at each receiver
f = (f + abs(min(f))*2);
f = f./(2*max(f));
% f = [0.5000    0.2693    0.2641    0.2602    0.2574    0.2555    0.2543    0.2538    0.2539    0.2546    0.2561    0.2583    0.2614    0.2656    0.2715    0.2795    0.2904    0.3048    0.3233    0.3458    0.3717].';
figure();
plot(theta,f,'k--')
title("Desired Beampattern")

N_a = length(xm);                       % only count remaining agents
N_s = length(rho);                      % want to sample transmitted strength to each receiver
mu = 2;                                 % 

%%
d = zeros(N_a,N_s);     % init fading parameters
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho(i,:).*[cos(theta(i,:));sin(theta(i,:))]);   % fading parameter 1 is distance
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));     % phase offset due to distance
    end
end

%%
figure();
scatter(xm,ym,'LineWidth',5) % plot workspace with agents
hold on
grid on
set(gca, 'LineWidth', 5, 'FontSize', 35)
scatter(rho(1).*cos(theta(1)),rho(1).*sin(theta(1)),'LineWidth',5)
scatter(rho(2:end).*cos(theta(2:end)),rho(2:end).*sin(theta(2:end)),'LineWidth',5)
xlabel('position x')
ylabel('position y')
axis('square')
xlim([-200 200])
ylim([-30 250])
legend('T_x','R_x','Adversary')