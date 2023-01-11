clc
clear all
rng('default');
c = @cmu.colors;

%%
gamma = 1;
mu = 2;
lambda = 3*10^8/(40*10^6);
k = 2*pi/lambda;
xm = lambda/2*[-2 -1 0 1 2]';
ym = zeros(5,1);
theta = [pi/6 pi/4]';
rho = 3;
f = [0.1 10]';

d = zeros(5,2);
zeta = zeros(5,2);
for m=1:5
    for i=1:2
        d(m,i) = norm([xm(m,:);ym(m,:)]-rho*[cos(theta(i,:));sin(theta(i,:))]);
        zeta(m,i) = k*(xm(m,:)*cos(theta(i,:)) + ym(m,:)*sin(theta(i,:)) + d(m,i));
    end
end

%%
T = 2*10^4;
a = zeros(5,T);
alpha = zeros(5,T);
grad = zeros(10,T-1);
grad_norm = zeros(1,T-1);
a(:,1) = rand(5,1);
alpha(:,1) = rand(5,1);
epsilon = 7*10^-3;

%%
for t=1:T-1
%     if t>1000
%         epsilon = 1/sqrt(t);
%     end
    u = zeros(5,2);
    v = zeros(5,2);
    for m=1:5
        for i=1:2
            u(m,i) = gamma/(d(m,i)^(mu/2))*cos(alpha(m,t)+zeta(m,i));
            v(m,i) = gamma/(d(m,i)^(mu/2))*sin(alpha(m,t)+zeta(m,i));
        end
    end
    ga = zeros(5,2);
    galpha = zeros(5,2);
    for i=1:2
        den = 0;
        for m=1:5
            den = den + a(m,t)*(u(m,i)+1i*v(m,i));
        end
        den = abs(den);
        num = den - f(i,:);
        ga(:,i) = num/den*((a(:,t)'*u(:,i))*u(:,i) + (a(:,t)'*v(:,i))*v(:,i));
        galpha(:,i) = num/den*(-(a(:,t)'*u(:,i))*a(:,t).*v(:,i) + (a(:,t)'*v(:,i))*a(:,t).*u(:,i));
    end
    a(:,t+1) = a(:,t) - epsilon*(ga(:,1)+ga(:,2));
    alpha(:,t+1) = alpha(:,t) - epsilon*(galpha(:,1)+galpha(:,2));
    grad(:,t) = [ga(:,1)+ga(:,2);galpha(:,1)+galpha(:,2)];
    grad_norm(:,t) = norm(grad(:,t));
    norm(grad(:,t))
end

%%
grad_norm(end)

figure(1)
plot(grad_norm, 'color', c('bright green'), 'LineWidth', 5)
xlabel('iteration t')
ylabel('||g(t)||')
set(gca, 'LineWidth', 5, 'FontSize', 35)
grid on
hold on
