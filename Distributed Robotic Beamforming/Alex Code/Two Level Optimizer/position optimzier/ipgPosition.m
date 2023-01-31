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
xOld = zeros(Na,T);
y = r(2,:).';
yOld = zeros(Na,T);
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
    u = zeros(Na,Ns);
    v = zeros(Na,Ns);
    dx = zeros(Na,Ns);
    dxx = zeros(Na,Ns);
    dy = zeros(Na,Ns);
    dyy = zeros(Na,Ns);
    dxy = zeros(Na,Ns);
    Lx = zeros(Na,Ns);
    Ly = zeros(Na,Ns);
    Lxx = zeros(Na);
    Lyy = zeros(Na);
    Lxy = zeros(Na);
    exponent = zeros(Na,Ns);
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

    % build gradient for amplitude and alpha (the fraction in eqns 10/11)
    % den is the rec'd AF at each rec w/o fading
    for rec = 1:Ns
        den = 0;
        for agent = 1:Na
            den = den + a(agent,t)*(u(agent,rec) + 1i*v(agent,rec));
            % derivative of Tx-Rx distance
            exponent(agent,rec) = exp(1j*(alpha(agent,t)* + k*(x(agent,:)*cos(tht(rec,:)) + y(agent,:)*sin(tht(rec,:)) + d(agent,rec))));
            dx(agent,rec) = (x(agent,:)-rho(rec,:)*cos(tht(rec,:)))*1/sqrt((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2);
            dy(agent,rec) = (y(agent,:)-rho(rec,:)*sin(tht(rec,:)))*1/sqrt((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2);
            dxx(agent,rec) = (1/sqrt((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2)) - (x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2*((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2)^(-3/2);
            dyy(agent,rec) = (1/sqrt((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2)) - (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2*((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2)^(-3/2);
            dxy(agent,rec) = (x(agent,:) - rho(rec,:)*cos(tht(rec,:)))*(y(agent,:) - rho(rec,:)*sin(tht(rec,:)))*((x(agent,:) - rho(rec,:)*cos(tht(rec,:)))^2 + (y(agent,:) - rho(rec,:)*sin(tht(rec,:)))^2)^(-3/2);
            % gradient of cost function
            % should be of dim'ality Na*Ns, summed along the dir of size Ns
            % to get gradient at each agent
            Lx(agent,rec) = (w(rec,:)/2)*((a(agent,t)*gamma(agent,rec)/d(agent,rec))*(k*(cos(tht(rec,:)) + dx(agent,rec))) - a(agent,t)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec)^(-2))*exponent(agent,rec);
            Ly(agent,rec) = (w(rec,:)/2)*((a(agent,t)*gamma(agent,rec)/d(agent,rec))*(k*(sin(tht(rec,:)) + dy(agent,rec))) - a(agent,t)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec)^(-2))*exponent(agent,rec);
            % hessian of cost function
            % should be of size Na*Na, used in multiplication - will need
            % to consist of two nested loops through each agent summed
            % across all Rxs
            % each of these 2nd derivative matricies should be of dim'ality
            % Na*Na, theyll stack into the full thing
            %             Lxx(agent,rec) = w(rec,:)*(-(a(agent,t)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec)^(-2))*(k*(cos(tht(rec,:)) + dx(agent,rec)))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*k*dxx(agent,rec))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*(k*(cos(tht(rec,:))+dx(agent,rec)))^2)*exponent(agent,rec) ...
            %                 -(a(agent,t)*gamma(agent,rec)*dxx(agent,rec)*d(agent,rec)^(-2))*exponent(agent,rec) ...
            %                 +2*(a(agent,t)*gamma(agent,rec)*dx(agent,rec)^2*d(agent,rec)^(-3))*exponent(agent,rec) ...
            %                 -(a(agent,rec)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec)^(-2))*k*(cos(tht(rec,:)) + dx(agent,rec))*exponent(agent,rec));
            %             Lyy(agent,rec) = w(rec,:)*(-(a(agent,t)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec)^(-2))*(k*(sin(tht(rec,:)) + dy(agent,rec)))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*k*dyy(agent,rec))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*(k*(sin(tht(rec,:))+dy(agent,rec)))^2)*exponent(agent,rec) ...
            %                 -(a(agent,t)*gamma(agent,rec)*dyy(agent,rec)*d(agent,rec)^(-2))*exponent(agent,rec) ...
            %                 +2*(a(agent,t)*gamma(agent,rec)*dy(agent,rec)^2*d(agent,rec)^(-3))*exponent(agent,rec) ...
            %                 -(a(agent,rec)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec)^(-2))*k*(sin(tht(rec,:)) + dy(agent,rec))*exponent(agent,rec));
            %             Lxy(agent,rec) = w(rec,:)*(-(a(agent,t)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec)^(-2))*(k*(sin(tht(rec,:)) + dy(agent,rec)))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*k*dxy(agent,rec))*exponent(agent,rec) ...
            %                 +(a(agent,t)*gamma(agent,rec)*d(agent,rec)^(-1)*(k*(sin(tht(rec,:))+dy(agent,rec)))*(k*(cos(tht(rec,:))+dx(agent,rec)))^2)*exponent(agent,rec) ...
            %                 -(a(agent,t)*gamma(agent,rec)*dxy(agent,rec)*d(agent,rec)^(-2))*exponent(agent,rec) ...
            %                 +2*(a(agent,t)*gamma(agent,rec)*dx(agent,rec)*dy(agent,rec)*d(agent,rec)^(-3))*exponent(agent,rec) ...
            %                 -(a(agent,rec)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec)^(-2))*k*(cos(tht(rec,:))+dx(agent,rec))*exponent(agent,rec));
            Lxx(1:Na,agent) = Lxx(1:Na,agent) + w(rec,:)*(-(a(:,t)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec).^(-2))*(k*(cos(tht(rec,:)) + dx(agent,rec)))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*k*dxx(agent,rec))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*(k*(cos(tht(rec,:))+dx(agent,rec))).^2)*exponent(agent,rec) ...
                -(a(:,t)*gamma(agent,rec)*dxx(agent,rec)*d(agent,rec).^(-2))*exponent(agent,rec) ...
                +2*(a(:,t)*gamma(agent,rec)*dx(agent,rec).^2*d(agent,rec).^(-3))*exponent(agent,rec) ...
                -(a(:,rec)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec).^(-2))*k*(cos(tht(rec,:)) + dx(agent,rec))*exponent(agent,rec));
            Lyy(1:Na,agent) = Lyy(1:Na,agent) + w(rec,:)*(-(a(:,t)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec).^(-2))*(k*(sin(tht(rec,:)) + dy(agent,rec)))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*k*dyy(agent,rec))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*(k*(sin(tht(rec,:))+dy(agent,rec))).^2)*exponent(agent,rec) ...
                -(a(:,t)*gamma(agent,rec)*dyy(agent,rec)*d(agent,rec).^(-2))*exponent(agent,rec) ...
                +2*(a(:,t)*gamma(agent,rec)*dy(agent,rec).^2*d(agent,rec).^(-3))*exponent(agent,rec) ...
                -(a(:,rec)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec).^(-2))*k*(sin(tht(rec,:)) + dy(agent,rec))*exponent(agent,rec));
            Lxy(1:Na,agent) = Lxy(1:Na,agent) + w(rec,:)*(-(a(:,t)*gamma(agent,rec)*dx(agent,rec)*d(agent,rec).^(-2))*(k*(sin(tht(rec,:)) + dy(agent,rec)))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*k*dxy(agent,rec))*exponent(agent,rec) ...
                +(a(:,t)*gamma(agent,rec)*d(agent,rec).^(-1)*(k*(sin(tht(rec,:))+dy(agent,rec)))*(k*(cos(tht(rec,:))+dx(agent,rec))).^2)*exponent(agent,rec) ...
                -(a(:,t)*gamma(agent,rec)*dxy(agent,rec)*d(agent,rec).^(-2))*exponent(agent,rec) ...
                +2*(a(:,t)*gamma(agent,rec)*dx(agent,rec)*dy(agent,rec)*d(agent,rec).^(-3))*exponent(agent,rec) ...
                -(a(:,rec)*gamma(agent,rec)*dy(agent,rec)*d(agent,rec).^(-2))*k*(cos(tht(rec,:))+dx(agent,rec))*exponent(agent,rec));
        end
        den = abs(den);
        AF(rec,t) = den;
        
        % numerator is error (diff btwn des'd and rec'd AF)
        % now calculate the gradient at each reciever

        % gx(:,rec) = w(rec,:)*0;
        % gy(:,rec) = w(rec,:)*(num/den*(-(a(:,t)'*u(:,rec))*a(:,t).*v(:,rec) + (a(:,t)'*v(:,rec))*a(:,t).*u(:,rec)));
    end
    Lx = abs(Lx); Ly = abs(Ly); Lxx = abs(Lxx); Lyy = abs(Lyy); Lxy = abs(Lxy);
    AFDB(:,t) = 20*log10(AF(:,t)/max(AF(:,t)));
    errorDB(t) = norm(AFDB(:,t) - 20*log10(f/max(f)),1)/Ns;


    %build hessian
    h(1:Na,1:Na) = Lxx;
    h(Na+1:end,1:Na) = Lxy.';
    h(1:Na,Na+1:end) = Lxy;
    h(Na+1:end,Na+1:end) = Lyy;
    % for agent = 1:Na
    %     for rec = 1:Ns
    %         den = 0;
    %         % rerun the denominator calculation on each agent
    %         for agent2 = 1:Na
    %             den = den + a(agent2,t)*(u(agent2,rec) + 1i*v(agent2,rec));
    %         end
    %         den = abs(den);
    %         num = den - f(rec,:);
    %
    %         % build hessian
    %         h(1:Na,agent) = h(1:Na,agent) + w(rec,:)*(num/den*(u(agent,rec)*u(:,rec)+v(agent,rec)*v(:,rec)) + f(rec,:)/den^3*(u(agent,rec)*u(:,rec)'*a(:,t)+v(agent,rec)*v(:,rec)'*a(:,t))*(u(:,rec)*u(:,rec)'*a(:,t)+v(:,rec)*v(:,rec)'*a(:,t)));
    %         h(1:Na,agent+Na) = h(1:Na,agent+Na) + w(rec,:)*(num/den*(-u(agent,rec)*a(:,t).*v(:,rec)-a(:,t)'*u(:,rec)*v(:,rec).*I(:,agent)+v(agent,rec)*a(:,t).*u(:,rec)+a(:,t)'*v(:,rec)*u(:,rec).*I(:,agent))...
    %             + f(rec,:)/den^3*((u(agent,rec)*u(:,rec)'*a(:,t)+v(agent,rec)*v(:,rec)'*a(:,t)))*(-a(:,t)'*u(:,rec)*a(:,t).*v(:,rec)+a(:,t)'*v(:,rec)*a(:,t).*u(:,rec)));
    %         h(agent+Na,1:Na) = h(1:Na,agent+Na)';
    %         h(Na+1:2*Na,agent+Na) = h(Na+1:2*Na,agent+Na) + w(rec,:)*(num/den*(a(agent,t)*v(agent,rec)*a(:,t).*v(:,rec)+a(agent,t)*u(agent,rec)*a(:,t).*u(:,rec)-a(:,t)'*u(:,rec)*a(:,t).*u(:,rec).*I(:,agent)-a(:,t)'*v(:,rec)*a(:,t).*v(:,rec).*I(:,agent))...
    %             + f(rec,:)/den^3*(-a(agent,t)*v(agent,rec)*a(:,t)'*u(:,rec)+a(agent,t)*u(agent,rec)*a(:,t)'*v(:,rec))*(-a(:,t)'*u(:,rec)*a(:,t).*v(:,rec)+a(:,t)'*v(:,rec)*a(:,t).*u(:,rec)));
    %     end
    % end


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
