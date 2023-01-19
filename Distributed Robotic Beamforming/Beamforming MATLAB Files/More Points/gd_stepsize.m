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
grad = zeros(2*N_a,T-1);
grad_norm = zeros(1,T-1);
a(:,1) = x(1:N_a,end);
alpha(:,1) = x(N_a+1:2*N_a,end);
epsilon = 10^-6;
beta = 1;

%%
for t=1:T-1
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
        ga(:,i) = num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i));
        galpha(:,i) = num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i));
    end
    h = zeros(2*N_a);
    for m=1:N_a
        for i=1:N_s
            den = 0;
            for j=1:N_a
                den = den + a(j,t)*(u(j,i)+1i*v(j,i));
            end
            den = abs(den);
            num = den - f(i,:);
            h(1:N_a,m) = h(1:N_a,m) + num/den*(u(m,i)*u(:,i)+v(m,i)*v(:,i)) + f(i,:)/den^3*(u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t))*(u(:,i)*u(:,i)'*a(:,t)+v(:,i)*v(:,i)'*a(:,t));
            h(1:N_a,m+N_a) = h(1:N_a,m+N_a) + num/den*(-u(m,i)*a(:,t).*v(:,i)-a(:,t)'*u(:,i)*v(:,i).*I(:,m)+v(m,i)*a(:,t).*u(:,i)+a(:,t)'*v(:,i)*u(:,i).*I(:,m))...
                + f(i,:)/den^3*((u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t)))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i));
            h(m+N_a,1:N_a) = h(1:N_a,m+N_a)';
            h(N_a+1:2*N_a,m+N_a) = h(N_a+1:2*N_a,m+N_a) + num/den*(a(m,t)*v(m,i)*a(:,t).*v(:,i)+a(m,t)*u(m,i)*a(:,t).*u(:,i)-a(:,t)'*u(:,i)*a(:,t).*u(:,i).*I(:,m)-a(:,t)'*v(:,i)*a(:,t).*v(:,i).*I(:,m))...
                + f(i,:)/den^3*(-a(m,t)*v(m,i)*a(:,t)'*u(:,i)+a(m,t)*u(m,i)*a(:,t)'*v(:,i))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i));  
        end
    end
    epsilon = 1/(max(eig(h))+beta); % change GD stepsize function of eigevals of hessian - optimal soln
    a(:,t+1) = a(:,t) - epsilon*sum(ga,2);
    alpha(:,t+1) = alpha(:,t) - epsilon*sum(galpha,2);
    grad(:,t) = [sum(ga,2);sum(galpha,2)];
    grad_norm(:,t) = norm(grad(:,t));
    disp(['iter: ',num2str(t),' error: ', num2str(norm(grad(:,t)))])
end

%%
figure(1)
plot(grad_norm, 'color', c('bright green'), 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on
