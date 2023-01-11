clc
c = @cmu.colors;

%%
d = zeros(N_a,N_s);
zeta = zeros(N_a,N_s);

%%
T = 1*10^4;
T_out = 1*10^2;
I = eye(N_a);
a = zeros(N_a,T);
alpha = zeros(N_a,T);
xm = zeros(N_a,T_out);
ym = zeros(N_a,T_out);
sep = lambda/4;
for i=1:N_a
    xm(i,1) = (mod(i,8)-4)*sep;
    ym(i,1) = (ceil(i/8)-5)*sep;
end
vmx = zeros(N_a,1);
vmy = zeros(N_a,1);
a(:,1) = 10*ones(N_a,1);
alpha(:,1) = 10*ones(N_a,1);
grad_norm = zeros(1,T-1);
grad_r_norm = zeros(1,T_out-1);
epsilon = 1*10^-4;
epsilon2 = 1*10^-5;
S = 1*I;
beta = 0;

%%
for outer=1:T_out
    for m=1:N_a
    for i=1:N_s
        d(m,i) = norm([xm(m,outer);ym(m,outer)]-rho*[cos(theta(i,:));sin(theta(i,:))]);
        zeta(m,i) = k*(xm(m,outer)*cos(theta(i,:)) + ym(m,outer)*sin(theta(i,:)) + d(m,i));
    end
    end
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
%     h = zeros(2*N_a);
%     for m=1:N_a
%         for i=1:N_s
%             den = 0;
%             for j=1:N_a
%                 den = den + a(j,t)*(u(j,i)+1i*v(j,i));
%             end
%             den = abs(den);
%             num = den - f(i,:);
%             h(1:N_a,m) = h(1:N_a,m) + w(i,:)*(num/den*(u(m,i)*u(:,i)+v(m,i)*v(:,i)) + f(i,:)/den^3*(u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t))*(u(:,i)*u(:,i)'*a(:,t)+v(:,i)*v(:,i)'*a(:,t)));
%             h(1:N_a,m+N_a) = h(1:N_a,m+N_a) + w(i,:)*(num/den*(-u(m,i)*a(:,t).*v(:,i)-a(:,t)'*u(:,i)*v(:,i).*I(:,m)+v(m,i)*a(:,t).*u(:,i)+a(:,t)'*v(:,i)*u(:,i).*I(:,m))...
%                 + f(i,:)/den^3*((u(m,i)*u(:,i)'*a(:,t)+v(m,i)*v(:,i)'*a(:,t)))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));
%             h(m+N_a,1:N_a) = h(1:N_a,m+N_a)';
%             h(N_a+1:2*N_a,m+N_a) = h(N_a+1:2*N_a,m+N_a) + w(i,:)*(num/den*(a(m,t)*v(m,i)*a(:,t).*v(:,i)+a(m,t)*u(m,i)*a(:,t).*u(:,i)-a(:,t)'*u(:,i)*a(:,t).*u(:,i).*I(:,m)-a(:,t)'*v(:,i)*a(:,t).*v(:,i).*I(:,m))...
%                 + f(i,:)/den^3*(-a(m,t)*v(m,i)*a(:,t)'*u(:,i)+a(m,t)*u(m,i)*a(:,t)'*v(:,i))*(-a(:,t)'*u(:,i)*a(:,t).*v(:,i)+a(:,t)'*v(:,i)*a(:,t).*u(:,i)));  
%         end
%     end
%     epsilon = 1/(max(eig(h))+beta);
    a(:,t+1) = a(:,t) - epsilon*sum(ga,2);
    alpha(:,t+1) = alpha(:,t) - epsilon*sum(galpha,2);
    grad(:,t) = [sum(ga,2);sum(galpha,2)];
    grad_norm(:,t) = norm(grad(:,t));
%     disp(['iter: ',num2str(t),' error: ', num2str(norm(grad(:,t)))])
   end
   gx = zeros(N_a,N_s);
   gy = zeros(N_a,N_s); 
   for m=1:N_a
        for i=1:N_s
            den = 0;
            for j=1:N_a
                den = den + a(j,t)*(u(j,i)+1i*v(j,i));
            end
            den = abs(den);
            num = den - f(i,:);
            k1 = -mu*a(m,t)*u(m,i)*(xm(m,outer)-rho*cos(theta(i,:))/(2*d(m,i)^2)) - a(m,t)*v(m,i)*k*(cos(theta(i,:))+(xm(m,outer)-rho*cos(theta(i,:)))/d(m,i));
            k2 = -mu*a(m,t)*v(m,i)*(xm(m,outer)-rho*cos(theta(i,:))/(2*d(m,i)^2)) + a(m,t)*u(m,i)*k*(cos(theta(i,:))+(xm(m,outer)-rho*cos(theta(i,:)))/d(m,i));
            k3 = -mu*a(m,t)*u(m,i)*(ym(m,outer)-rho*sin(theta(i,:))/(2*d(m,i)^2)) - a(m,t)*v(m,i)*k*(sin(theta(i,:))+(ym(m,outer)-rho*sin(theta(i,:)))/d(m,i));
            k4 = -mu*a(m,t)*v(m,i)*(ym(m,outer)-rho*sin(theta(i,:))/(2*d(m,i)^2)) + a(m,t)*u(m,i)*k*(sin(theta(i,:))+(ym(m,outer)-rho*sin(theta(i,:)))/d(m,i));
            gx(m,i) = num/den*(a(:,t)'*u(:,i)*k1 + a(:,t)'*v(:,i)*k2);
            gy(m,i) = num/den*(a(:,t)'*u(:,i)*k3 + a(:,t)'*v(:,i)*k4);
        end
   end
   xm(:,outer+1) = xm(:,outer) - epsilon2*(sum(gx,2)-vmx);
   ym(:,outer+1) = ym(:,outer) - epsilon2*(sum(gy,2)-vmy);
   vmx = vmx - epsilon2*S*(xm(:,outer+1)-xm(:,1));
   vmy = vmy - epsilon2*S*(ym(:,outer+1)-ym(:,1));
   grad_r(:,outer) = [sum(gx,2)-vmx;sum(gy,2)-vmy];
   grad_r_norm(:,outer) = norm(grad_r(:,outer));
   disp(['out_iter: ',num2str(outer),' error: ', num2str(norm(grad_r(:,outer)))])
%    disp(['out_iter: ',num2str(outer),' x1: ', num2str(xm(1,outer))])
   
end

%%
% figure(1)
% plot(grad_norm, 'color', c('bright green'), 'LineWidth', 5)
% xlabel('iteration t')
% ylabel('||g(t)||')
% set(gca, 'LineWidth', 5, 'FontSize', 35)
% grid on
% hold on
