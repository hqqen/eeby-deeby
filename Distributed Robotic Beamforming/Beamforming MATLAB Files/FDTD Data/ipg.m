clc
rng('default')

%%
T = 2500;
I = eye(N_a);
a = zeros(N_a,T);
alpha = zeros(N_a,T);
grad = zeros(2*N_a,T-1);
grad_norm = zeros(1,T-1);
af =  zeros(N_s,1);
af_db = zeros(N_s,T-1);
a(:,1) = 10*rand(N_a,1);
alpha(:,1) = 0*rand(N_a,1);
epsilon1 = 1*10^-5;
epsilon2 = 1;
K = 0*ones(2*N_a);
beta = 0.001;
gamma = ones(N_a,N_s);
u_ch = zeros(N_a,N_s);
v_ch = zeros(N_a,N_s);

w = 1./sqrt(f);
w = w/sum(w);

%%
for t=1:T-1
    u = zeros(N_a,N_s);
    v = zeros(N_a,N_s);
    for m=1:N_a
        for i=1:N_s
            u(m,i) = gamma(m,i)/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i));
            v(m,i) = gamma(m,i)/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
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
%         af(i,:) = den;
        num = den - f(i,:);
        ga(:,i) = w(i,:)*(num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i)));
        galpha(:,i) = w(i,:)*(num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i))); 
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
            h(1:N_a,m) = h(1:N_a,m) + w(i,:)*(num/den*(u(m,i)*u(:,i)+v(m,i)*v(:,i)) + f(i,:)/den^3*(u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t))*(u(:,i)*u(:,i)'*a(:,t)+v(:,i)*v(:,i)'*a(:,t)));
            h(1:N_a,m+N_a) = h(1:N_a,m+N_a) + w(i,:)*(num/den*(-u(m,i)*a(:,t).*v(:,i)-a(:,t)'*u(:,i)*v(:,i).*I(:,m)+v(m,i)*a(:,t).*u(:,i)+a(:,t)'*v(:,i)*u(:,i).*I(:,m))...
                + f(i,:)/den^3*((u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t)))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));
            h(m+N_a,1:N_a) = h(1:N_a,m+N_a)';
            h(N_a+1:2*N_a,m+N_a) = h(N_a+1:2*N_a,m+N_a) + w(i,:)*(num/den*(a(m,t)*v(m,i)*a(:,t).*v(:,i)+a(m,t)*u(m,i)*a(:,t).*u(:,i)-a(:,t)'*u(:,i)*a(:,t).*u(:,i).*I(:,m)-a(:,t)'*v(:,i)*a(:,t).*v(:,i).*I(:,m))...
                + f(i,:)/den^3*(-a(m,t)*v(m,i)*a(:,t)'*u(:,i)+a(m,t)*u(m,i)*a(:,t)'*v(:,i))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));  
        end
    end
    epsilon1 = 1/(max(eig(h))+beta);
    K = K - epsilon1*(h*K+beta*K-eye(2*N_a));
    grad(:,t) = [sum(ga,2);sum(galpha,2)];
    x = [a(:,t);alpha(:,t)] - epsilon2*K*grad(:,t);
    a(:,t+1) = x(1:N_a,:);
    alpha(:,t+1) = x(N_a+1:2*N_a,:);
    grad_norm(:,t) = norm(grad(:,t));
    
    for m=1:N_a
        for i=1:N_s
            u_ch(m,i) = gamma_ch(m,i)*cos(alpha(m,t)+zeta(m,i));
            v_ch(m,i) = gamma_ch(m,i)*sin(alpha(m,t)+zeta(m,i));
%             u_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i));
%             v_ch(m,i) = gamma_ch(m,i)/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
        end
    end
    for i=1:N_s
        den_ch = 0;
        for m=1:N_a
            den_ch = den_ch + a(m,t)*(u_ch(m,i)+1i*v_ch(m,i));
        end
        den_ch = abs(den_ch);
        af(i,:) = den_ch;
    end
    af_db(:,t) = 20*log10(af/max(af));
    error_fdb = norm(20*log10(af/max(af))-20*log10(f/max(f)),1)/N_s;
    error_f = norm(af-f,1)/N_s;
    disp(['iter: ',num2str(t),' error_db: ',num2str(error_fdb),' error: ',num2str(error_f)])
    
%     if error_fdb < 1.8
%         beta = 0.5;
%         epsilon2 = 1;
%     end
   
end

%%
% figure(1)
% plot(theta,af_db(:,1), 'g', 'LineWidth', 5)
% xlabel('theta (radian)')
% ylabel('radiation pattern (dB)')
% set(gca, 'LineWidth', 5, 'FontSize', 35)
% grid on
% hold on

%%
for t=1:t-1
error_fdb(:,t) = norm(af_db(:,t)-20*log10(f/max(f)),1)/N_s;
end
[p,q]=min(error_fdb);
% plot(theta,af_db(:,q), 'g', 'LineWidth', 3)
% hold on
% plot(theta,20*log10(f/max(f)), '--k', 'LineWidth', 3)