clear all; close all;
% algorithm params
f0 = 40e6;           % beam freq
Na = 8;             % num Tx
Nb = 9;              % num Tx for building signal
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

xm = [1 2 3 4 5 6 7 8].';
ym = [1 2 3 4 5 6 7 8].';
% xm = [15    15    15    15    15    15   105   105   105   105   105   105   195   195   195   196   195   195   285   285   285   285   285   285   375   375   375   378   375   375   465   465   465   465   465   465].';
% ym = [15   105   195   285   375   465    16   105   195   286   375   465    15   105   200   285   371   465    15   105   195   285   375   465    16   105   195   285   375   465    15   105   195   285   376   465].';
 r = [xm, ym].';

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

%% get and plot AF
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

%% build Rx Posns and weights
rho = [];
for i = 1:max(size(theta))

    rho(i) = 5*d;
    w(i) = 1/sqrt(f(i));

end
rho = rho.'; w = w.';
w = w/sum(w);

%% build initla agent weights and amplitudes
a = []; alpha = [];
for i = 1:Na    
    a(i) = b(i);
    alpha(i) = b(i)*pi/8;
end
%% plot Tx/Rx setup
figure(2);
scatter(r(1,:), r(2,:), 'b', 'LineWidth', 3); hold on;
scatter(rho.*cos(theta),rho.*sin(theta), 'k', 'LineWidth', 3)
legend("Tx", "Rx")

%% get cost
L = getCost(Na,Ns,xm,ym,theta,w,rho,alpha,a,f)
