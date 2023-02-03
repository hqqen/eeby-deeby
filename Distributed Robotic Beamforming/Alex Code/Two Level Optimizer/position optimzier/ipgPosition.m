function [amp, alpha, AF, noisyAF] = ipgPhaseMag(T,f,r,rho,tht,a0,alpha0,w,f0)
% run IPGD-BF on the phase and magntiude of a set of transmitters
% inputs:
% T - number of optimizer steps to take
% f - desired beampattern at each sampled point
% r - 2xn vector showing the position of n agents
% rho - reciever distance
% tht - angular position of each reciever (Rx_i @ (rho,tht)_i in polar)
% a0 - initial cplx amplitude guess
% alpha0 - initial phase offsest guess
% f0 - beam frequency (probably 40e6)
% w - vector of error weights (1/rx)
% returns:
% amp - 1xNa vector of cplx transmitter amplitudes
% alpha - 1xNa vector of transmitte rphase offsets

% setup algorithm parameters
eps1 = 1e-2; % GD rate for precond matrix
eps2 = 1;    % parameter GD rate
beta = .1;   % stabilization parameter

% build other algorithm parameters
Na = length(r);     % num agents
Ns = length(rho);   % num points to sample over
gamma = ones(Na,Ns);% placeholder no noise gamma
mu = 2;             % path loss exponent
a = ones(Na,T).*reshape(a0,[Na,1]);
alpha = ones(Na,T).*reshape(alpha0,[Na,1]);
r = reshape(r,2,max(size(r)));
x = r(1,:).';         % split agent posn vector r into x and y posns
x0 = x;
y = r(2,:).';
y0 = y;
lambda = 3e8/f0;
k = (2*pi)/lambda;
I = eye(Na);
K = 0*ones(2*Na);
tht = reshape(tht, max(size(tht)), 1);
w = reshape(w,max(size(w)), 1);
f = reshape(f, max(size(f)), 1);


% initialize channel fading parameters
d = zeros(Na,Ns);       % transmitter - reciever distance
zeta = zeros(Na,Ns);    % exponent for finding AF
for agent = 1:Na
    for rec = 1:Ns
        d(agent,rec) = norm( [x(agent,:);y(agent,:)] - rho(rec,:).*[cos(tht(rec,:));sin(tht(rec,:))]);
        zeta(agent,rec) = k*x(agent,:)*cos(tht(rec,:)) + k*y(agent,:)*sin(tht(rec,:)) + k*d(agent,rec);
    end
end

% optimizer loop
for t = 1:T
    %clear parameters udated stepwise
    % u = zeros(Na,Ns);
    % v = zeros(Na,Ns);
    % dx = zeros(Na,Ns);
    % dxx = zeros(Na,Ns);
    % dy = zeros(Na,Ns);
    % dyy = zeros(Na,Ns);
    % dxy = zeros(Na,Ns);
    % Lx = zeros(Na,Ns);
    % Ly = zeros(Na,Ns);
    % Lxx = zeros(Na);
    % Lyy = zeros(Na);
    % Lxy = zeros(Na);
    % exponent = zeros(Na,Ns);
%     gx = zeros(Na,Ns);
%     gy = zeros(Na,Ns);
    h = zeros(2*Na);

    % update the noiseless channel parameters to send to server
    for agent = 1:Na
        for rec = 1:Ns
            u(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*cos(alpha(agent,t) + zeta(agent,rec));
            v(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*sin(alpha(agent,t) + zeta(agent,rec));
        end
    end

    % start building gradient and hessian for algorithm
    % get 1st deriv of dist to Rxs and dist traveled by agents
    % update distance traveled by each agent
    p = zeros(Na,1);
    for agent = 1:Na
        p(agent) = sqrt((x(agent) - x0(agent))^2 + (y(agent) - y0(agent))^2);
    end

    % build exponent
    for agent = 1:Na
        for rec = 1:Ns
            sigma(agent,rec) = alpha(agent) + k*(x(agent)*cos(tht(rec)) + y(agent)*sin(rec) + d(agent,rec));
        end
    end

    dx = zeros(Na,1); dy = zeros(Na,1); px = zeros(Na,1); py = zeros(Na,1); sigmax = zeros(Na); sigmay = zeros(Na);
    for rec = 1:Ns
        for agent = 1:Na
            dx(agent) = dx(agent) + (x(agent) - rho(rec)*cos(tht(rec)))/d(agent,rec);
            dy(agent) = dy(agent) + (y(agent) - rho(rec)*sin(tht(rec)))/d(agent,rec);
            px(agent) = (x(agent) - x0(agent))/p(agent);
            py(agent) = (y(agent) - y0(agent))/p(agent);
            sigmax(agent,rec) = k*cos(tht(rec)) + dx(agent);
            sigmay(agent,rec) = k*sin(tht(rec)) + dy(agent);
        end
    end
    % get 2nd deriv of dist to Rxs and dist traveled by agents
    % because the distance from agent agent to a Rx and its dist traveled
    % is not effected by the motion of other agents these 2nd derives
    % end up being Na*Na diag'l matrices
    dxx = eye(Na); dyy = eye(Na); dxy = eye(Na); pxx = eye(Na); pyy = eye(Na); pxy = eye(Na); sigmaxx = zeros(Na); sigmayy = zeros(Na); sigmaxy = zeros(Na);
    for rec = 1:Ns
        for agent = 1:Na
            dxx(agent,agent) = dxx(agent,agent) + (1/d(agent,rec)) - (x(agent) - rho(rec)*cos(tht(rec)))*dx(agent)/(d(agent,rec)^2);
            dyy(agent,agent) = dyy(agent,agent) + (1/d(agent,rec)) - (y(agent) - rho(rec)*sin(tht(rec)))*dy(agent)/(d(agent,rec)^2);
            dxy(agent,agent) = dxy(agent,agent) - dx(agent)*dy(agent)/d(agent,rec);
            pxx(agent,agent) = pxx(agent,agent) + (1/p(agent)) - (x(agent) - x0(agent))*px(agent)/(p(agent)^2);
            pyy(agent,agent) = pyy(agent,agent) + (1/p(agent)) - (y(agent) - y0(agent))*py(agent)/(p(agent)^2);
            pxy(agent,agent) = pxy(agent,agent) - px(agent)*py(agent)/p(agent);
            % are these three right? they might need to be diagonal
            % sigmaxx(agent,rec) = dxx(agent,agent);
            % sigmayy(agent,rec) = dyy(agent,agent);
            % sigmaxy(agent,rec) = dxy(agent,agent);
        end
    end
    simgaxx = dxx; sigmayy = dyy; sigmaxy = dxy;
    % build first derivatives for grad calculation
    Lx = zeros(Na,Ns); Ly = zeros(Na,Ns);
    for rec = 1:Ns
        for agent = 1:Na
            Lx(agent) = Lx(agent) + (a(agent)*gamma(agent,rec)/d(agent,rec))*exp(1j*sigma(agent,rec))*(-dx(agent)/d(agent,rec) + sigmax(agent,rec));
            Ly(agent) = Ly(agent) + (a(agent)*gamma(agent,rec)/d(agent,rec))*exp(1j*sigma(agent,rec))*(-dy(agent)/d(agent,rec) + sigmay(agent,rec));
        end
    end
    % build second derivs for hessian calculaiton
    for rec = 1:Ns
        for agent = 1:Na
            % Lxx(agent,agent) = 
            % Lyy(agnet,agent) = 
            % Lxy(agent,agent) = 
        end
    end

    % update GD parameter(s)
    eps1 = 1/(max(eig(h)) + beta);
    % run GD update
    K = K - eps1*(h*K+beta*K-eye(2*Na));         % mass update preconditioner (using eqn 5, not 12)
    grad(:,t) = [sum(Lx,2);sum(Ly,2)];            % agent-wise gradient is summed
    g = [x(:);y(:)] - eps2*K*grad(:,t);   % gradient update (eqn 4)
    x(:) = g(1:Na);                            % split gradient update to amplitude and phase
    y(:) = g(Na+1:end);
    grad_norm(:,t) = norm(grad(:,t));                 % take l2 norm of grad

    % now calculate the true rec'd AF under the effect of noise
    noisyGamma = normrnd(1,.1,Na,Ns);
    for agent = 1:Na
        for rec = 1:Ns
            noisyU(agent,rec) = noisyGamma(agent,rec)/(d(agent,rec)^(mu/2))*cos(alpha(agent,t)+zeta(agent,rec));
            noisyV(agent,rec) = noisyGamma(agent,rec)/(d(agent,rec)^(mu/2))*sin(alpha(agent,t)+zeta(agent,rec));
        end
    end
    % find true array factor
    for rec = 1:Ns
        den = 0;
        for agent = 1:Na
            den = den + a(agent,t)*(noisyU(agent,rec)+1i*noisyV(agent,rec));
        end
        noisyAF(rec,t) = abs(den);

    end
    noisyAFDB(:,t) = 20*log10(noisyAF(:,t)/max(noisyAF(:,t)));
    noisyErrorDB(t) = norm(noisyAFDB(:,t)-20*log10(f/max(f)),1)/Ns;
    if isnan(noisyErrorDB(t))
        sprintf("Algorithm has Diverged!")
        break
    end
    if ~mod(t,100)

        figure(20202);
        plot(tht, noisyAFDB(:,t), 'g', 'LineWidth', 3); hold on
        plot(tht, 20*log10(f/max(f)), '--k', 'LineWidth', 3);
        xlabel("\theta (rad)")
        ylabel("Beampattern (dB)")
        title('Noisy Beampattern')
        grid on;
        set(gca, 'FontSize', 12);
        hold off;

        figure(30303);
        plot(movmean(noisyErrorDB(1:t),25), 'LineWidth', 3)
        title("Noisy Windowed Average of Error (k = 25)")
        ylabel("||E(t)|| (dB)")
        xlabel("t")
        grid on;
        set(gca, 'FontSize', 12)

        % figure(40404);
        % plot(tht,AFDB(:,t), 'g', 'LineWidth', 3); hold on;
        % plot(tht, 20*log10(f/max(f)), '--k', 'LineWidth', 3);
        % xlabel("\theta (rad)")
        % ylabel("Beampattern (dB)")
        % title('Noiseless Beampattern')
        % grid on;
        % % set(gca, 'FontSize', 24);
        % hold off;
        %
        % figure(50505);
        % plot(movmean(errorDB(1:t),25), 'LineWidth', 3)
        % title("Noiseless Windowed Average of Error (k = 25)")
        % ylabel("||E(t)|| (dB)")
        % xlabel("t")
        % grid on;
        % set(gca, 'FontSize', 24)

        figure(60606); hold off
        scatter(x,y,'b','LineWidth',4); hold on
        scatter(rho.*cos(tht),rho.*sin(tht),'k','LineWidth',4)
        title("Agent Positions")
        xlabel("x (m)"); ylabel("y (m)");
%         exportgraphics(gcf,'agentMovementNoSBL.gif','Append',true);
        grid on


    end
end
figure(6060606)                                   % plot
plot(tht,noisyAFDB(:,1), 'g', 'LineWidth', 5); hold on
plot(tht,noisyAFDB(:,T/10), 'r', 'LineWidth', 5)
plot(tht,noisyAFDB(:,T/5), 'b', 'LineWidth', 5)
plot(tht,noisyAFDB(:,T/2), 'c', 'LineWidth', 5)
plot(tht,noisyAFDB(:,T), 'm', 'LineWidth', 5)
set(gca, 'FontSize', 35)
plot(tht,20*log10(f/max(f)), '--k', 'LineWidth', 3)
title("Beampattern Matching with Noise")
legend("t = 0", "t = " + string(T/10), "t = " + string(T/5), "t = " + string(T/2), "t = " + string(T), "Desired", "Location", "Best")
xlabel('\theta (radian)')
ylabel('radiation pattern (dB)')
set(gca, 'LineWidth', 5, 'FontSize', 12)
drawnow
grid on

drawnow
end
