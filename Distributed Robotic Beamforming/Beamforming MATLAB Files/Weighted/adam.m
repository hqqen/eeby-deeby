clc
c = @cmu.colors;

%%
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);
for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho*[cos(theta(i,:));sin(theta(i,:))]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end

%%
T = 5*10^4;
I = eye(N_a);
a = zeros(N_a,T);
alpha = zeros(N_a,T);
a(:,1) = 10*ones(N_a,1);
alpha(:,1) = 10*ones(N_a,1);
grad_norm = zeros(1,T-1);
epsilon = 0.005;
beta = 0;
mean1 = 0;
mean2 = 0;
var1 = 0;
var2 = 0;

%%
for t=1:T-1
%     if t>10 && grad_norm(:,t-1)>grad_norm(:,t-2)
%         epsilon1 = 0.0001;
%     end
    u = zeros(N_a,N_s);
    v = zeros(N_a,N_s);
    for m=1:N_a
        for i=1:N_s
            u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i));
            v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
        end
    end
    ga = zeros(N_a,N_s);
    galpha = zeros(N_a,N_s);
    for i=1:N_s
        den = 0;
        for m=1:N_a
            den = den + a(m,t)*(u(m,i)+1i*v(m,i));
        end
        den = abs(den);
        num = den - f(i,:);
        ga(:,i) = w(i,:)*(num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i)));
        galpha(:,i) = w(i,:)*(num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i)));
    end
    epsilon = 1/sqrt(t);
    if t>0.5*10^4
        epsilon = 0.1/sqrt(t);
    elseif t>1.5*10^4
        epsilon = 0.05/sqrt(t);
    elseif t>2.5*10^4
        epsilon = 0.01/sqrt(t);
    end
    mean1 = 0.9*mean1 + 0.1*sum(ga,2);
    mean2 = 0.9*mean2 + 0.1*sum(galpha,2);
    var1 = 0.999*var1 + 0.001*sum(ga,2).^2;
    var2 = 0.999*var2 + 0.001*sum(galpha,2).^2;
    a(:,t+1) = a(:,t) - epsilon*mean1./(sqrt(var1)+10^-8)*sqrt(1-0.999^t)/(1-0.9^t);
    alpha(:,t+1) = alpha(:,t) - epsilon*mean2./(sqrt(var2)+10^-8)*sqrt(1-0.999^t)/(1-0.9^t);
    grad(:,t) = [sum(ga,2);sum(galpha,2)];
    grad_norm(:,t) = norm(grad(:,t));
    disp(['iter: ',num2str(t),' error: ', num2str(norm(grad(:,t)))])
end

%%
figure(1)
plot(grad_norm, 'color', 'm', 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on
