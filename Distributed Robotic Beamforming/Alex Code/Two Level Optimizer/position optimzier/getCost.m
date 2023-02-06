function [L] = getCost(Na,Ns,x,y,tht,w,rho,alpha,a,f)
lambda = 3e8/4e7;
k = (2*pi)/lambda;
mu = 2;
x = reshape(x,Na,1);
y = reshape(y,Na,1);
tht = reshape(tht,Ns,1);
alpha = reshape(alpha,Na,1);
a = reshape(a,Na,1);
L = 0;

d = zeros(Na,Ns);       % transmitter - reciever distance
zeta = zeros(Na,Ns);    % exponent for finding AF
for agent = 1:Na
    for rec = 1:Ns
        d(agent,rec) = norm( [x(agent,:);y(agent,:)] - rho(rec,:).*[cos(tht(rec,:));sin(tht(rec,:))]);
        zeta(agent,rec) = k*x(agent,:)*cos(tht(rec,:)) + k*y(agent,:)*sin(tht(rec,:)) + k*d(agent,rec);
    end
end

% optimizer loop

%clear parameters udated stepwise
u = zeros(Na,Ns);
v = zeros(Na,Ns);
ga = zeros(Na,Ns);
galpha = zeros(Na,Ns);
h = zeros(2*Na);

% update the noiseless channel parameters to send to server
gamma = ones(Na,Ns);
for agent = 1:Na
    for rec = 1:Ns
        u(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*cos(alpha(agent,1) + zeta(agent,rec));
        v(agent,rec) = gamma(agent,rec)./(d(agent,rec).^(mu/2)).*sin(alpha(agent,1) + zeta(agent,rec));
    end
end

% build gradient for amplitude and alpha (the fraction in eqns 10/11)
% den is the rec'd AF at each rec w/o fading
for rec = 1:Ns
    den = 0;
    for agent = 1:Na
        den = den + a(agent,1)*(u(agent,rec) + 1i*v(agent,rec));
    end
    den = abs(den);
    AF(rec,1) = den;
    % numerator is error (diff btwn des'd and rec'd AF)
    num = den - f(rec, :);
    % now calculate the gradient at each reciever
    ga(:,rec) = w(rec,:)*(num/den*((a(:,1)'*u(:,rec))*u(:,rec) + (a(:,1)'*v(:,rec))*v(:,rec)));
    galpha(:,rec) = w(rec,:)*(num/den*(-(a(:,1)'*u(:,rec))*a(:,1).*v(:,rec) + (a(:,1)'*v(:,rec))*a(:,1).*u(:,rec)));
    L = L + (w(rec)/2)*abs(f(rec) - abs(AF(rec,1)));
end
% AFDB(:,1) = 20*log10(AF(:,1)/max(AF(:,1)));
% errorDB(1) = norm(AFDB(:,1) - 20*log10(f/max(f)),1)/Ns;

end