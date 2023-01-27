% Parameters of desired beam
f = 40e6;                                     % Frequency in Hz
N_d = 21;                                       % Number of elements of array
lambda = 3e8/f;                               % wave length
d = .5*lambda;                                 % Inter-element spacing
k = (2*pi)./lambda;                           % wave number
% channel parameters
gamma = 1;
mu = 2;
nu = 1;
rho = 1.5*d;

s = 360;
%rho = [50*lambda/2*ones(1,1);50*lambda/2*ones(s/2,1);60*lambda/2*ones(s/2,1)].'; % get distance from transmitter to receiver
theta = pi./15;                         % angle of allied reciever
theta_sec = linspace(0,6*pi./3,s);  % divide workspace into 30deg sectors
theta = [theta_sec'].';             % make list of all receiver posns

% s = 20;                             % numer of sectors to look at
% f = [5;1*ones(s,1)]/10;             % set magnitude of des'd beampower at each receiver
% % rho = [50*lambda/2*ones(1,1);50*lambda/2*ones(s/2,1);60*lambda/2*ones(s/2,1)]; % get distance from transmitter to receiver
% theta = pi/15;                         % angle of allied reciever
% theta_sec = linspace(pi./3,3*pi./3,s);  % divide workspace into 30deg sectors
% theta = [theta;theta_sec'].';             % make list of all receiver posns

rm_d = []; Am_d = []; alpha_d = [];
c = 1:floor(N_d/2); c = [c, c(end)+1, flip(c)];
beam = [];
figure();
for s = 0:.01:2*pi
    for i = 1:N_d

        rm_d(:,i) = [0;d * c(i)];
        Am_d(i) = c(i);
        alpha_d(i) = s*c(i);%c(i)*pi/8;%pi/(180/N_d) * (i - ceil(N_d/2));

    end

    % Array factor (AFd) of desired beam
    AFd = GetArrayFactor(N_d, rm_d, Am_d, alpha_d, gamma, mu, nu, rho, theta); % call provided code to find AF
    AF = abs(AFd);
    h = polar(theta,AF,'r'); % plot
    set(h,'LineWidth',2);
    box on
    % hold on
    drawnow

    beam(end+1) = sum(AFd./max(AFd));
end

figure();
plot(0:.01:2*pi, beam);
title("Beam Directivity")

AF = abs(AFd);
f = AF;
f = (f + abs(min(f))*2);
f = f./(2*max(f)).';