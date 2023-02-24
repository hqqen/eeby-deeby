clear all;
Na = 6; Ns = 4;
%% setup symvars
syms x y a alpha x0 y0 [Na 1] real
syms rho tht [Ns 1] real
syms rMin eps1 eps2 real
mu = 2; k = (2*pi)/(3e8/4e7);
sympref('FloatingPointOutput', true); % this should go in the file invoking this, not here!
digits(16)
%% reshape inputs
Ns = numel(rho);
% Na = numel(xTrue);
rho = reshape(rho,[numel(rho),1]);
tht = reshape(tht,[numel(tht),1]);
a = reshape(a,[numel(a),1]);
alpha = reshape(alpha,[numel(alpha),1]);
% xTrue = reshape(xTrue,[numel(xTrue),1]);
% yTrue = reshape(yTrue,[numel(yTrue),1]);
% syms rho tht [4 1] real
% syms k mu real positive

%% initialize storage
d = sym(zeros(Na,Ns));
zeta = sym(zeros(Na,Ns));
u = sym(zeros(Na,Ns));
v = sym(zeros(Na,Ns));
dr = sym(zeros(2*Na,Ns));
ys = sym(zeros(Ns,1));
dAgent = sym(zeros(Na,Na));
p = sym(zeros(Na,Na));

%% build array of distance based penalties 
% need to build in travel dist penalties and inter agent spacing penalties
r = [x,y]; r0 = [x0, y0];
for i = 1:Na
    for j = 1:Na
        % dAgent(i,j) = eps1*norm(r(i,:) - r(j,:))^2; % travel dist
        pAgent(i,j) = eps2*log(norm(r(i,:) - r(j,:)) - rMin); % interagent spacing
        pAgent(i,i) = eps1*norm(r(i,:) - r0(i,:)); % travel distance
    end
end
pAgent = sum(pAgent,2); % sum total of added cost at each agent

%% build terms for the first derivative
r = [x;y];
for i = 1:max(size(rho))
    for m = 1:max(size(x))
        d(m,i) = sqrt((x(m) - rho(i)*cos(tht(i)))^2 + (y(m) - rho(i)*sin(tht(i)))^2);
        zeta(m,i) = k*(x(m)*cos(tht(i)) + y(m)*sin(tht(i)) + d(m,i));
    end
    u(:,i) = (1./d(:,i)).*cos(alpha + zeta(:,i));
    v(:,i) = (1./d(:,i)).*sin(alpha + zeta(:,i));
    % find apparent AF at rec i
    ys(i) = sqrt(u(:,i).'*a*a.'*u(:,i) + v(:,i).'*a*a.'*v(:,i));
    % find the gradient w/r.t. position, sum across the 1st dim so it's
    % [2Na*Ns] (i.e. it shows the change in cost at each Rx for moving in a
    % given dimension)
    dr(:,i) = sum((1/ys(i))*(a.'*u(:,i)*jacobian(u(:,i),r) + a.'*v(:,i)*jacobian(v(:,i),r)).',2);
end

%% sum cost for each coordiante (i.e get a 2Na*1 array)
% now add the distance traveled penalties and interagent spacing
dp = sum(jacobian(pAgent,r).',2);
% dp = sign(dp).*min(abs(dp),10);
dr = sum(dr,2) + dp;
%% get second derivative
% ddr = subs(jacobian(sum(dr,2),r),[x,y],[xTrue,yTrue]);
fprintf("Finding Hessian...\n")
tic
ddr = matlabFunction(vpa(jacobian(dr,r)),'Optimize',false,'Vars',{x,y,x0,y0,a,alpha,rho,tht,rMin,eps1,eps2},'file','getHessianFullCost');
matlabFunction(vpa(dp),'Vars',{r,r0,rMin,eps1,eps2},'file','getCostDerivs','Optimize',false)
% ddr = double(subs(jacobian(sum(dr,2),r),[x,y],[xTrue,yTrue]));
fprintf("Hessian Found! \n")
toc

