clear all; close all;
% algorithm params
f0 = 40e6;           % beam freq
Na = 33;             % num Tx
Nb = 36;              % num Tx for building signal
lambda = 3e8/f0;     % wavelength
d = lambda/2;       % interagent spacing

% channel params (dont change!)
gamma = 1;
mu = 2;
nu = 1;
rhoCh = 3*d/2;

%% build array factor and Tx parameters (one call is plotting, one is to get feasible profile)
NsCh = 360; % sample 360 points to build complete AF profile
Ns = 30;    % sample 20 ponts for the algorithm to optimize over
thetaCh = linspace(0,2*pi,NsCh);
theta = linspace(0,7*pi/8,Ns);
r = []; a0 = []; alpha0 = [];
c = 1:ceil(Na/2); c = [c, flip(c)];
% setup disance bounds for Rxs
dMin = -10*d; dMax = 10*d;
for i = 1:Na

    r(:,i) = (dMax - dMin)*rand(2,1) + dMin;
    % a0(i) = c(i)^3;
    % alpha0(i) = c(i)*pi/8;

end
% make intial amplitude and phase guess
a0 = 100*rand(Na,1);
alpha0 = rand(Na,1);
% get desired array factor
b = 1:ceil(Nb/2); b = [b, flip(b)];
for i = 1:Nb

    rb(:,i) = [0, i*d/2];
    a(i) = b(i);
    alpha(i) = b(i)*pi/8;

end
AFd = GetArrayFactor(Nb,rb,a,alpha,gamma,mu,nu,rhoCh,theta);
f = abs(AFd).';
% sample the entire array factor
AFdCh = GetArrayFactor(Nb,rb,a,alpha,gamma,mu,nu,rhoCh,thetaCh);
fCh = abs(AFdCh).';
%plot array factor
figure(1);
plot(thetaCh,20*log10(fCh)/max(fCh),'r', 'LineWidth', 5); hold on;
plot(theta,20*log10(f)/max(f),'g', 'LineWidth', 5);
set(gca,'FontSize',12)
grid on;
legend("Full AF", "Sampled AF")
xlabel("\theta (rad)")

%% get distance to each Rx and set error weights
rho = [];
for i = 1:max(size(theta))

    rho(i) = 50*d;
    w(i) = 1/sqrt(f(i));

end
rho = rho.'; w = w.';
w = w/sum(w);


%% plot Tx/Rx setup
figure(2);
scatter(r(1,:), r(2,:), 'b', 'LineWidth', 3); hold on;
scatter(rho.*cos(theta),rho.*sin(theta), 'k', 'LineWidth', 3)
legend("Tx", "Rx")

%% run optimizer
% first run SBL to prune agents
KK = Main_SBL(r, f.', Na, theta);

% plot the agent array after pruning
figure(41);
scatter(r(1,abs(KK)>0),r(2,abs(KK)>0),'b','LineWidth',4); hold on;
scatter(r(1,abs(KK)==0),r(2,abs(KK)==0),'r','LineWidth',4)
title("Agents After Pruning")
legend("Kept Agents", "Pruned Agents")
grid on

% rebuild agent array to only have agents which survived pruning
rPruned = r(:, abs(KK) > 0);
a0 = abs(KK(abs(KK) > 0));
alpha0 = angle(KK(abs(KK) > 0));

% then run IPG to further optimize weights
ipgPhaseMag(1e4,f,rPruned,rho,theta,a0,alpha0,w,f0)