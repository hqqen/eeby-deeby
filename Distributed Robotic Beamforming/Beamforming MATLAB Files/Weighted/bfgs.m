clc

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
gnew = zeros(2*N_a,T-1);
grad_norm = zeros(1,T-1);
a(:,1) = 10*ones(N_a,1);
alpha(:,1) = 10*ones(N_a,1);
% alpha = 10^-3;
Binv = eye(2*N_a);

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
        ga(:,i) = w(i,:)*(num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i)));
        galpha(:,i) = w(i,:)*(num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i))); 
    end
    grad(:,t) = [sum(ga,2);sum(galpha,2)];
    s = -Binv*grad(:,t);
%     delta = backtrack(@F,@DF,x(:,t),s,alpha);
    delta = 10^-6;
    x = [a(:,t);alpha(:,t)] + delta*s;
    a(:,t+1) = x(1:N_a,:);
    alpha(:,t+1) = x(N_a+1:2*N_a,:);
    
    u = zeros(N_a,N_s);
    v = zeros(N_a,N_s);
    for m=1:N_a
        for i=1:N_s
            u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t+1)+zeta(m,i));
            v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t+1)+zeta(m,i));
        end
    end
    ga = zeros(N_a,N_s);
    galpha = zeros(N_a,N_s);
    for i=1:N_s
        den = 0;
        for m=1:N_a
            den = den + a(m,t+1)*(u(m,i)+1i*v(m,i));
        end
        den = abs(den);
        num = den - f(i,:);
        ga(:,i) = w(i,:)*(num/den*((a(:,t+1)'*u(:,i))*u(:,i) + (a(:,t+1)'*v(:,i))*v(:,i)));
        galpha(:,i) = w(i,:)*(num/den*(-(a(:,t+1)'*u(:,i))*a(:,t+1).*v(:,i) + (a(:,t+1)'*v(:,i))*a(:,t+1).*u(:,i))); 
    end
    gnew(:,t) = [sum(ga,2);sum(galpha,2)];
    
    y = gnew(:,t) - grad(:,t);
    Binv = Binv + (s'*y+y'*Binv*y)*(s*s')/(s'*y)^2 -(Binv*y*s'+s*y'*Binv)...
        /(s'*y);
    grad_norm(:,t) = norm(gnew(:,t));
    disp(['iter: ',num2str(t),' error: ', num2str(norm(gnew(:,t)))])
%     if t>100 && grad_norm(:,t)>10^4
%         break;
%     end
end

%%
figure(1)
plot(grad_norm(1:915), 'color', 'm', 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on