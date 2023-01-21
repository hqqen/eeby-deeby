clc
c = @cmu.colors;

%%
d = zeros(N_a,N_s); % initialize channel fading parameters parameters d & zeta
zeta = zeros(N_a,N_s);
for m=1:N_a                                                                             % for each agent...
    for i=1:N_s                                                                         % ...at each sampled angle
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho*[cos(theta(i,:));sin(theta(i,:))]);         % d is the distance from agent m to receiver i
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));     % channel fading (which is unknown to agents) is modeled as a function of wave number and distance
    end
end

%%
T = 1.5*10^4;                   % run for 1.5e4 iterations
I = eye(N_a);                   % create identity matrix of size nAgents x nAgents
a = zeros(N_a,T);               % init agent amplitudes
alpha = zeros(N_a,T);           % init agent phases
grad = zeros(2*N_a,T-1);        % init gradient storage (2*nAgents, need to store phase & magn gradient)
grad_norm = zeros(1,T-1);       % init storage for l2 norm of gradient
a(:,1) = 10*ones(N_a,1);        % init amplitudes at 10 for each agent
alpha(:,1) = 10*ones(N_a,1);    % init phases at 10 for each agent
epsilon1 = 1*10^-5;             % preconditioner gradient descent rate - epsilon in paper
epsilon2 = 1;                   % parameter gradient descent rate
K = 0*ones(2*N_a);              % init preconditioner K
beta = 0;                       % stabilization parameter
vFade = zeros(N_a, N_s);        % init storage for u and v with fading taken into account
uFade = zeros(N_a, N_s);

%%
for t=1:T-1
    u = zeros(N_a,N_s);         % init vectors to be updated agent-wise
    v = zeros(N_a,N_s);
    for m=1:N_a
        for i=1:N_s
            u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i)); % for each agent run the real valued update equation
            v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i)); % ...and the imaginary valued update equation
        end
    end
    ga = zeros(N_a,N_s);        % init gradient storage for amplitude
    galpha = zeros(N_a,N_s);    % ...and for phase
    for i=1:N_s
        den = 0;                                    % init update eqn denominator
        for m=1:N_a
            den = den + a(m,t)*(u(m,i)+1i*v(m,i));  % for each agent add the rec'd beampattern params (in the absence of fading?)
        end
        den = abs(den);                                                                                 % get magnitude of rec'd beam
        num = den - f(i,:);                                                                             % numerator is difference in rec'd and des'd beam
        ga(:,i) = num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i));                          % eqn 10 (w/o a(t))
        galpha(:,i) = num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i));     % eqn 11 (w/o alpha(t))
    end
    h = zeros(2*N_a);                                   % init hessian
    for m=1:N_a                                         % for each agent...
        for i=1:N_s                                     % ...at each sampled angle
            den = 0;                                    % init denominator
            for j=1:N_a                                 % for each agent at each sampled angle...
                den = den + a(j,t)*(u(j,i)+1i*v(j,i));  % ...sum all rec'd beams
            end
            den = abs(den);                                                                                                                                                                 % find magnitude of rec'd signal
            num = den - f(i,:);                                                                                                                                                             % numerator is difference in rec'd and des'd power
            % compute hessian (done by amplitude then by phase, then by the
            % mirrored terms then finally by cross terms)
            h(1:N_a,m) = h(1:N_a,m) + num/den*(u(m,i)*u(:,i)+v(m,i)*v(:,i)) + f(i,:)/den^3*(u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t))*(u(:,i)*u(:,i)'*a(:,t)+v(:,i)*v(:,i)'*a(:,t));
            h(1:N_a,m+N_a) = h(1:N_a,m+N_a) + num/den*(-u(m,i)*a(:,t).*v(:,i)-a(:,t)'*u(:,i)*v(:,i).*I(:,m)+v(m,i)*a(:,t).*u(:,i)+a(:,t)'*v(:,i)*u(:,i).*I(:,m))...
                + f(i,:)/den^3*((u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t)))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i));
            h(m+N_a,1:N_a) = h(1:N_a,m+N_a)';
            h(N_a+1:2*N_a,m+N_a) = h(N_a+1:2*N_a,m+N_a) + num/den*(a(m,t)*v(m,i)*a(:,t).*v(:,i)+a(m,t)*u(m,i)*a(:,t).*u(:,i)-a(:,t)'*u(:,i)*a(:,t).*u(:,i).*I(:,m)-a(:,t)'*v(:,i)*a(:,t).*v(:,i).*I(:,m))...
                + f(i,:)/den^3*(-a(m,t)*v(m,i)*a(:,t)'*u(:,i)+a(m,t)*u(m,i)*a(:,t)'*v(:,i))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i));
        end
    end
    %     epsilon1 = 0.5/(max(eig(h))+beta);
    epsilon1 = 10^-7;  % update GD weights
    epsilon2 = .8;     % "
    Kold = K;
    K = K - epsilon1*(h*K+beta*K-eye(2*N_a));       % preconditioner update (eqn 5)
    grad(:,t) = [sum(ga,2);sum(galpha,2)];          % get gradient for current timestep
    x = [a(:,t);alpha(:,t)] - epsilon2*K*grad(:,t); % gradient update equation (eqns 10-11)
    a(:,t+1) = x(1:N_a,:);                          % update amplitudes
    alpha(:,t+1) = x(N_a+1:2*N_a,:);                % update phases
    grad_norm(:,t) = norm(grad(:,t));               % l2 norm of grad for plot
    disp(['iter: ',num2str(t),' error: ', num2str(norm(grad(:,t)))])

    %plot rec'd AF
    if mod(t,100) == 0
        gamma_ch = normrnd(1,0.1,N_a,N_s);                % channel fading parameter is randomized at each step
        for m=1:N_a                                                                  % for each agent
            for i=1:N_s                                                              % at each sampled angle
                u_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i)); % agent-wise update vector w/fading
                v_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
            end
        end
        for i=1:N_s                                                 % at each sampled angle
            den_ch = 0;                                             % reset rec'd AF (?)
            for m=1:N_a                                             % for each agent
                den_ch = den_ch + a(m,t)*(u_ch(m,i)+1i*v_ch(m,i));  % update rec'd AF (eqn 9)
            end
            den_ch = abs(den_ch);                                   % magnitude of AF
            af(i,t) = den_ch;                                       % store rec'd(?) AF
        end
        af_db(:,t) = 20*log10(af(:,t)/max(af(:,t)));                % convert norm'd AF to dB
        error_fdb = norm(af_db(:,t)-20*log10(f/max(f)),1)/N_s;      % find error btwn norm'd rec'd AF and norm'd des'd AF

        figure(2)                                   % plot
        plot(theta,af_db(:,t), 'g', 'LineWidth', 5)
        hold on
        plot(theta,20*log10(f/max(f)), '--k', 'LineWidth', 3)
        hold off
        xlabel('theta (radian)')
        ylabel('radiation pattern (dB)')
        set(gca, 'LineWidth', 5, 'FontSize', 35)
        drawnow
        grid on
    end

    % if t>1700 && grad_norm(:,t)>grad_norm(:,t-1)
    %     save('x_ipg.mat','x');
    %     save('K.mat','Kold');
    %     break;
    % end
end

%%
figure(1) % plot l2 norm of grad
plot(grad_norm, 'color', 'blue', 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on
