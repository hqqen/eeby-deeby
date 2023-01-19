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
T = 10^4;                       % run for 1e4 iterations
a = zeros(N_a,T);               % intialize agent amplitudes
alpha = zeros(N_a,T);           % intialize agent phases offsets
grad = zeros(2*N_a,T-1);        % intialize gradient storage (dim 2*nAgents b/c each grad is magn; phase) 
grad_norm = zeros(1,T-1);       % intialize storage for l2 norm of gradient
a(:,1) = 10*ones(N_a,1);        % for the first GD pass every agent amplitude is set to 10
alpha(:,1) = 10*ones(N_a,1);    % ...and every agent phase offset is also set to 10
epsilon = 10^-6;                % gradient descent rate

%%
for t=1:T-1
%     if t>2000
%         epsilon = 10^-3;
%     end
    u = zeros(N_a,N_s);         % initialize the vectors to be computer agent-wise 
    v = zeros(N_a,N_s);
    for m=1:N_a
        for i=1:N_s
            u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i)); % for each agent use local info to compute the real part of its beamforming magnitude
            v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i)); % ...as well as the imaginary part of the beamforming magnitude
        end
    end
    ga = zeros(N_a,N_s);                                              % gradient descent storage for agent amplitude
    galpha = zeros(N_a,N_s);                                          % ...and for agent phase offsets
    for i=1:N_s
        den = 0;                                                      % initialize denominator
        for m=1:N_a
            den = den + a(m,t)*(u(m,i)+1i*v(m,i));                    % update denominator agent-wise by adding amplitude*AF (cplx)
        end
        den = abs(den);
        num = den - f(i,:);                                           % difference in predicted AF magnitude and actual AF magnitude (given in f)
        ga(:,i) = num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i));                      % get gradient for magnitude for current iteration
        galpha(:,i) = num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i)); % get gradient for phase for current iteration
    end
    a(:,t+1) = a(:,t) - epsilon*sum(ga,2);                            % update the amplitude for the next sim step
    alpha(:,t+1) = alpha(:,t) - epsilon*sum(galpha,2);                % update phase offset for the next step
    grad(:,t) = [sum(ga,2);sum(galpha,2)];                            % sum gradient agent-wise across all sampled angles w/in beam
    grad_norm(:,t) = norm(grad(:,t));                                 % take l2 norm of current timestep gradient
    disp(['iter: ',num2str(t),' error: ', num2str(norm(grad(:,t)))])
end

%%
figure(1)   % plot magnitude of gradient as function of iteration t
plot(grad_norm, 'color', c('bright green'), 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on
