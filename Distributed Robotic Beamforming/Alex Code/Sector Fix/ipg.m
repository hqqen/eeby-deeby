clc
% rng('default')

%%
T = 1e4;                   % run for 1e3 iterations
I = eye(N_a);               % for later use
a = zeros(N_a,T);           % array of amplitude guesses
alpha = zeros(N_a,T);       % array of phase guesses
grad = zeros(2*N_a,T-1);    % gradient storage
grad_norm = zeros(1,T-1);   % l2 norm gradient storage
af =  zeros(N_s,1);         % calc'd array factor
af_db = zeros(N_s,T-1);     % calc'd array factor in dB
a(:,1) = 100*rand(N_a,1);   % the first amplitude guess is randomly dist'd from 0-100
alpha(:,1) = rand(N_a,1);   % initial phase guess is randomly dist'd from 0-1
epsilon1 = 1e-2;            % preconditioner descent rate (eqn 5)
epsilon2 = 1;               % gradient descent rate (delta, eqn 4)
K = 0*ones(2*N_a);          % init preconditioner
beta = 0.1;                 % local preconditioner update rate (eqns 5, 12)
gamma = ones(N_a,N_s);      % init channel fading params
u_ch = zeros(N_a,N_s);      % init agentwise post-fading update vectors
v_ch = zeros(N_a,N_s);
error_fdb = zeros(1,T);     %error storage

% w = 1./sqrt(f);             % init penalty weights
% w = (w/sum(w));               % normalize penalty weights
 w = .01*ones(size(f));

%%
for t=1:T-1                % for each iteration
    u = zeros(N_a,N_s);    % init storage for agent-wise update vectors
    v = zeros(N_a,N_s);
    % step 1
    for m=1:N_a            % for each agent
        for i=1:N_s        % for each sampled angle
            u(m,i) = gamma(m,i)/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i)); % run agent-wise update equation (eqns 7, 8)
            v(m,i) = gamma(m,i)/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
        end
    end
    % step 2 happens between lines here (theres no need to actually send
    % information btwn agents)
    ga = zeros(N_a,N_s);        % gradient storage for amplitude
    galpha = zeros(N_a,N_s);    % gradient storage for phase
    %step 3
    for i=1:N_s                                     % for each sampled angle
        den = 0;                                    % reset denominator for
        for m=1:N_a                                 % for each agent
            den = den + a(m,t)*(u(m,i)+1i*v(m,i));  % sum across all agent instances of eqn 9
        end
        % end step 3 - at this point the server would push u(t), v(t),
        % y(t) = sum(a(i)*(u(i) + j*v(i))) and K(t) to agents
        % step 4 happens btwn lines here similar to step 2
        % begin step 5
        den = abs(den);        % take magnitude of AF
        % af(i,:) = den;         % store clac'd array factor
        num = den - f(i,:);    % numerator is difference in rec'd and des'd AF (eqn 10)
        % these two lines are the summands of eqns 10 & 11
        ga(:,i) = w(i,:)*(num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i)));                         % calcualte amplitude gradient (step 5, eqn 10)
        galpha(:,i) = w(i,:)*(num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i)));    % calculate phae gradient (step 5, eqn 11)
    end

    h = zeros(2*N_a);                                   % init hessian
    for m=1:N_a                                         % for each agent
        for i=1:N_s                                     % at each sampled angle
            den = 0;                                    % reset denominator
            for j=1:N_a                                 % for each agent
                den = den + a(j,t)*(u(j,i)+1i*v(j,i));  % rebuild the fraction for finding hessian of eqn 3
            end
            den = abs(den);                             % get magn of AF (for use in eqns 10, 11)
            num = den - f(i,:);                         % numerator is difference of rec'd and des'd AF (for use in eqn 10)
            % calculate hessian (mirrored terms first then nonmirrored terms)
            h(1:N_a,m) = h(1:N_a,m) + w(i,:)*(num/den*(u(m,i)*u(:,i)+v(m,i)*v(:,i)) + f(i,:)/den^3*(u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t))*(u(:,i)*u(:,i)'*a(:,t)+v(:,i)*v(:,i)'*a(:,t)));
            h(1:N_a,m+N_a) = h(1:N_a,m+N_a) + w(i,:)*(num/den*(-u(m,i)*a(:,t).*v(:,i)-a(:,t)'*u(:,i)*v(:,i).*I(:,m)+v(m,i)*a(:,t).*u(:,i)+a(:,t)'*v(:,i)*u(:,i).*I(:,m))...
                + f(i,:)/den^3*((u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t)))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));
            h(m+N_a,1:N_a) = h(1:N_a,m+N_a)';
            h(N_a+1:2*N_a,m+N_a) = h(N_a+1:2*N_a,m+N_a) + w(i,:)*(num/den*(a(m,t)*v(m,i)*a(:,t).*v(:,i)+a(m,t)*u(m,i)*a(:,t).*u(:,i)-a(:,t)'*u(:,i)*a(:,t).*u(:,i).*I(:,m)-a(:,t)'*v(:,i)*a(:,t).*v(:,i).*I(:,m))...
                + f(i,:)/den^3*(-a(m,t)*v(m,i)*a(:,t)'*u(:,i)+a(m,t)*u(m,i)*a(:,t)'*v(:,i))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));
        end
    end

    epsilon1 = 1/(max(eig(h))+beta);                % pick optimal descent rate for preconditioner
%     epsilon2 = .1/(max(eig(h)) - min(eig(h)));       % not in original code
    if t > 1000
        epsilon2 = 1/(ceil(t/1000));
        epsilon1 = .75/(ceil(t/1000)*(max(eig(h))+beta));
        beta = .1/ceil(t/1000);
    end
    K = K - epsilon1*(h*K+beta*K-eye(2*N_a));         % mass update preconditioner (using eqn 5, not 12)
    grad(:,t) = [sum(ga,2);sum(galpha,2)];            % agent-wise gradient is summed
    x = [a(:,t);alpha(:,t)] - epsilon2*K*grad(:,t);   % gradient update (eqn 4)
    a(:,t+1) = x(1:N_a,:);                            % split gradient update to amplitude and phase
    alpha(:,t+1) = x(N_a+1:2*N_a,:);
    grad_norm(:,t) = norm(grad(:,t));                 % take l2 norm of grad

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

    error_fdb(t) = norm(af_db(:,t)-20*log10(f/max(f)),1)/N_s;      % find error btwn norm'd rec'd AF and norm'd des'd AF in dB
    if isnan(error_fdb(t))
        break                                                      % if the algorithm ever returns an invalid error break
    end


    %     error_f = norm(af-f,1)/N_s;                               % find error not in dB
    disp(['iter: ',num2str(t),' error_db: ',num2str(error_fdb(t))])
    if ~mod(t,100)
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

        figure(3);
        plot(movmean(error_fdb(1:t),25))
        title("Windowed Average of Error in dB (k = 25)")
    end
end


figure(); hold on
for i = 1:707
    scatter(i,min(error_fdb(1:i)))
end

%%
% error = zeros(1,t-1);
% for i=1:t-1
%     error(:,i) = norm(af_db(:,i)-20*log10(f/max(f)),1)/N_s;
% end
% [p,q]=min(error);
% figure(3)
% plot(error)
%%
% figure(3)                                   % plot
% plot(theta,af_db(:,1), 'g', 'LineWidth', 5)
% xlabel('theta (radian)')
% ylabel('radiation pattern (dB)')
% set(gca, 'LineWidth', 5, 'FontSize', 35)
% grid on
% hold on
figure(4)
plot(movmean(error_fdb,100))
title("Windowed Average of Error in dB (k = 10)")
