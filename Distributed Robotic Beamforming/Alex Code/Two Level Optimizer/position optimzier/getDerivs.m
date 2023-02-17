function ddr = getDerivs(k,mu,rho,tht,Na,a,alpha,xTrue,yTrue)
%% setup symvars
syms x y [Na 1] real
sympref('FloatingPointOutput', true); % this should go in the file invoking this, not here!
%% reshape inputs
Ns = numel(rho);
% Na = numel(xTrue);
rho = reshape(rho,[numel(rho),1]);
tht = reshape(tht,[numel(tht),1]);
a = reshape(a,[numel(a),1]);
alpha = reshape(alpha,[numel(alpha),1]);
xTrue = reshape(xTrue,[numel(xTrue),1]);
yTrue = reshape(yTrue,[numel(yTrue),1]);
% syms rho tht [4 1] real
% syms k mu real positive

%% initialize storage
d = sym(zeros(Na,Ns));
zeta = sym(zeros(Na,Ns));
u = sym(zeros(Na,Ns));
v = sym(zeros(Na,Ns));
dr = sym(zeros(2*Na,Ns));
ys = sym(zeros(Ns,1));


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
%% get second derivative
% ddr = subs(jacobian(sum(dr,2),r),[x,y],[xTrue,yTrue]);
fprintf("Finding Hessian...\n")
tic
ddr = matlabFunction(jacobian(sum(dr,2),r),'Optimize',false,'Vars',{x,y},'file','getHessianOpt');
% ddr = double(subs(jacobian(sum(dr,2),r),[x,y],[xTrue,yTrue]));
fprintf("Hessian Found! \n")
toc
end

