clear all;
% test with 6 agents and 4 sampled points
Na = 4; Ns = 2;
% syms gamma [Na Ns] real
% syms d(x1, x2, x3, x4, y1, y2, y3, y4) sigma(x1, x2, x3, x4, y1, y2, y3, y4) real
syms f tht w rho [1 Ns] real
syms x y a alpha x0 y0 [1 Na] real
syms k real
L = 0; AF = sym(zeros(1,Ns));
gamma = ones(Na,Ns);
for agent = 1:Na
    % p(agent) = sqrt((x(agent) - x0(agent))^2 + (y(agent) - y0(agent))^2);
    % recAF = 0;
    for rec = 1:Ns
        d(agent,rec) = sqrt((x(agent) - rho(rec)*cos(tht(rec)))^2 + (y(agent) - rho(rec)*sin(tht(rec)))^2);
        % dx(agent) = dx(agent) + (x(agent) - rho(rec)*cos(tht(rec)))/d(agent,rec);
        % dy(agent) = dy(agent) + (y(agent) - rho(rec)*sin(tht(rec)))/d(agent,rec);
        % dxx(agent,)
        sigma(agent,rec) = alpha(agent) + k*(x(agent)*cos(tht(rec)) + y(agent)*sin(rec) + d(agent,rec));
        % recAF = recAF + a(agent)
        % L = L + (w(rec)/2) * abs(f(rec) - );%abs((a(agent)*gamma(agent,rec)/d(agent,rec))*exp(1j*(alpha(agent) + k*(x(agent)*cos(tht(rec)) + y(agent)*sin(tht(rec)) + d(agent,rec))))) + p(agent));
        % AF(rec) = AF(rec) + (a(agent)*gamma(agent,rec)/d(agent,rec))*exp(1j*(alpha(agent) + k*(x(agent)*cos(tht(rec)) + y(agent)*sin(tht(rec)) + d(agent,rec))));
    end
end

for rec = 1:Ns
    AFrec = 0;
    for agent = 1:Na
        % AFrec = AFrec + (a(agent)*gamma(agent,rec)/d(agent,rec))*exp(1j*(alpha(agent) + k*(x(agent)*cos(tht(rec)) + y(agent)*sin(tht(rec)) + d(agent,rec))));
        AFrec = AFrec + (a(agent)*gamma(agent,rec)/d(agent,rec))*exp(i*sigma(agent,rec));
    end
    L = L + (w(rec)/2)*(f(rec) - (AFrec));
end