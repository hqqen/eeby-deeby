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
a = reshape(a0,[Na,1]);
alpha = reshape(alpha0,[Na,1]);
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

    % initialize channel fading parameters
    d = zeros(Na,Ns);       % transmitter - reciever distance
    zeta = zeros(Na,Ns);    % exponent for finding AF
    for agent = 1:Na
        for rec = 1:Ns
            d(agent,rec) = norm( [x(agent,:);y(agent,:)] - rho(rec,:).*[cos(tht(rec,:));sin(tht(rec,:))]);
            zeta(agent,rec) = k*x(agent,:)*cos(tht(rec,:)) + k*y(agent,:)*sin(tht(rec,:)) + k*d(agent,rec);
        end
    end

    % update the noiseless channel parameters to send to server
    for agent = 1:Na
        for rec = 1:Ns
            u(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*cos(alpha(agent) + zeta(agent,rec));
            v(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*sin(alpha(agent) + zeta(agent,rec));
        end
    end

    % find derivatives of channel params
    % these vectors need to be 2*Na each (2 rows, Na cols) so their square
    % is Na*Na to keep dimensionality correct
    % the actual multiplicaiton will be done row-wise to give a 2Na*2Na
    % matrix (shaped like a hessian I think - [u^2, u*v; u*v, v^2]?)
    du = zeros(Na,2); dv = zeros(Na,2);
    for agent = 1:Na
        for rec = 1:Ns
            dux(agent,rec) = a(agent)*-(mu*cos(alpha(agent) + zeta(agent,rec))*(2*x(agent) - 2*rho(rec)*cos(tht(rec)))/(4*d(agent,rec)^(2*mu/2)) - (sin(alpha(agent) + zeta(agent,rec))*(k*cos(tht(rec)) + ((k*(2*x(agent) - 2*rho(rec)*cos(tht(rec))))/(d(agent,rec)^(mu/2)))))/d(agent,rec)^(mu/2));
            duy(agent,rec) = a(agent)*-(mu*cos(alpha(agent) + zeta(agent,rec))*(2*y(agent) - 2*rho(rec)*sin(tht(rec)))/(4*d(agent,rec)^(2*mu/2)) - (sin(alpha(agent) + zeta(agent,rec))*(k*sin(tht(rec)) + ((k*(2*y(agent) - 2*rho(rec)*sin(tht(rec))))/(d(agent,rec)^(mu/2)))))/d(agent,rec)^(mu/2));
            dvx(agent,rec) = a(agent)*-(mu*sin(alpha(agent) + zeta(agent,rec))*(2*x(agent) - 2*rho(rec)*cos(tht(rec)))/(4*d(agent,rec)^(2*mu/2)) + (cos(alpha(agent) + zeta(agent,rec))*(k*cos(tht(rec)) + ((k*(2*x(agent) - 2*rho(rec)*cos(tht(rec))))/(d(agent,rec)^(mu/2)))))/d(agent,rec)^(mu/2));
            dvy(agent,rec) = a(agent)*-(mu*sin(alpha(agent) + zeta(agent,rec))*(2*y(agent) - 2*rho(rec)*sin(tht(rec)))/(4*d(agent,rec)^(2*mu/2)) + (cos(alpha(agent) + zeta(agent,rec))*(k*sin(tht(rec)) + ((k*(2*y(agent) - 2*rho(rec)*sin(tht(rec))))/(d(agent,rec)^(mu/2)))))/d(agent,rec)^(mu/2));
        end
    end
    % build gradients
    for rec = 1:Ns
        % find the rec'd AF and the respective error
        recAF = 0;
        for agent = 1:Na
            recAF = recAF + a(agent)*(u(agent,rec) + 1i*v(agent,rec));
        end
        recAF = abs(recAF);
        recErr = recAF - f(rec);
        %  build the actual gradient
        gx(:,rec) = w(rec,:)*((recErr/recAF)*((a(:)'*u(:,rec))*dux(:,rec) + (a(:)'*v(:,rec))*dvx(:,rec)));
        gy(:,rec) = w(rec,:)*((recErr/recAF)*((a(:)'*u(:,rec))*duy(:,rec) + (a(:)'*v(:,rec))*dvy(:,rec)));
    end
    % get 2nd derivaives of u and v
    % for agent = 1:Na
    %     for rec = 1:Ns
    %         uxx(agent,rec) = -a(agent)*((-mu*sin(a(agent) + zeta(agent,rec))*(k*cos(tht(rec)) + k*(x(agent) - (rho(rec)*cos(tht(rec)))/d(agent,rec))*(2*x(agent) - 2*rho(agent)*cos(tht(agent))) + 2*mu*cos(alpha(agent) + zeta(agent,rec))))/(4*d(agent,rec)^(2+mu/2)) ...
    %                          - (mu*cos(alpha(agent) + zeta(agent,rec))*(2*x(agent)*2*rho(rec)*cos(tht(rec)))*(2+mu/2)*(d^(1+mu/2))*(x(agent)-rho(rec)*cos(tht(rec)))/d(agent,rec))/(4*d(agent,rec)^(2+mu/2))^2) ...
    %                          + (((cos(alpha(agent) + zeta(agent,rec))*(k*cos(tht(rec)) + k*(x(agent - rho(rec)*cos(tht(rec))))/d(agent,rec))*(k*cos(tht(rec)) + k*(x(agent - rho(rec)*cos(tht(rec))))/d(agent,rec)) + sin(alpha(agent) + zeta(agent,rec))*(4*k*d(agent,rec) - k*(2*x(agent) - 2*rho(rec)*cos(tht(rec)))^2/(2*d(agent,rec))^2))))
    % 


    % update GD parameter(s)
    eps1 = 1/(max(eig(h)) + beta);
    % run GD update
    % K = K - eps1*(h*K+beta*K-eye(2*Na));         % mass update preconditioner (using eqn 5, not 12)
    K = eye(2*Na);
    grad(:,t) = [sum(gx,2);sum(gy,2)];            % agent-wise gradient is summed
    g = [x(:);y(:)] - eps2*K*grad(:,t);   % gradient update (eqn 4)
    x(:) = g(1:Na);                            % split gradient update to amplitude and phase
    y(:) = g(Na+1:end);
    grad_norm(:,t) = norm(grad(:,t));                 % take l2 norm of grad

    % now calculate the true rec'd AF under the effect of noise
    noisyGamma = normrnd(1,.1,Na,Ns);
    for agent = 1:Na
        for rec = 1:Ns
            noisyU(agent,rec) = noisyGamma(agent,rec)/(d(agent,rec)^(mu/2))*cos(alpha(agent)+zeta(agent,rec));
            noisyV(agent,rec) = noisyGamma(agent,rec)/(d(agent,rec)^(mu/2))*sin(alpha(agent)+zeta(agent,rec));
        end
    end
    % find true array factor
    for rec = 1:Ns
        den = 0;
        for agent = 1:Na
            den = den + a(agent)*(noisyU(agent,rec)+1i*noisyV(agent,rec));
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
