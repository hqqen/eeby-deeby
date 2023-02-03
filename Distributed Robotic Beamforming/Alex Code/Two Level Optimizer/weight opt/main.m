clear all; close all;
% algorithm params
f0 = 40e6;           % beam freq
Na = 33;             % num Tx
Nb = 39;              % num Tx for building signal
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
dMin = -5*d; dMax = 5*d;
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

% fix agent positions for testing
% xm = [1 2 3 4 5 6 7 8];
% ym = [1 2 3 4 5 6 7 8];
% xm = [15    15    15    15    15    15   105   105   105   105   105   105   195   195   195   196   195   195   285   285   285   285   285   285   375   375   375   378   375   375   465   465   465   465   465   465].';
% ym = [15   105   195   285   375   465    16   105   195   286   375   465    15   105   200   285   371   465    15   105   195   285   375   465    16   105   195   285   375   465    15   105   195   285   376   465].';
% r = [17.8322933400355	-0.810765386231672	-21.3761523566949	-11.6788480302041	-33.6841594058367	2.56912937569614	14.7063260239211	28.7923920958415	18.9233175448750	0.0665033069661050	6.76950862273910	35.4567139328590	-21.3289268565745	14.3041254035042	-4.01563958904831	10.5907954188767	16.0332305933528	29.4094234368190	-3.87393110193165	-24.7019573296387	6.96004886600709	34.7303265472191	-6.75746466421664	-31.8122880099696	-34.2542097836698	-19.2888182674276	-9.45625498505702	-15.3429343333715	22.3970807913259	-28.5515868781852	-24.2562681236973	-21.2354683549748	-14.9421399861831;
% -4.10100805746026	17.2269070618586	-34.7056200528825	10.4731831386929	-8.35952040009143	3.23459835454643	-34.5838616527315	22.8667651292164	21.9622443988624	-5.82281781505257	24.4555767539464	-0.271499296048241	13.5434755776113	9.43213117074116	4.54604659411793	-26.8310368011747	-30.7192855667246	0.0383565053338373	19.8152667131031	-0.655249563279348	-36.2076274300987	14.6029018137280	-4.30717633051243	-22.8323125313450	-34.5576262344101	9.51360786498231	-16.4662104223187	8.47941766080180	21.1184827428183	-25.7174817537144	-23.1980085629934	10.1657697594529	-17.2528149355206];
% r = [xm; ym];
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
rPruned = r;r(:, abs(KK) > 0);
a0 = 100*rand(Na,1);abs(KK(abs(KK) > 0));
alpha0 = rand(Na,1);angle(KK(abs(KK) > 0));
% 
% % then run IPG to further optimize weights
ipgPhaseMag(1e4,f,rPruned,rho,theta,a0,alpha0,w,f0)