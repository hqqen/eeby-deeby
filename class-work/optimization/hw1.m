% homework 1
%% problem 2.3
fprintf("Q1: \n")
syms x [1 3] real
f = exp(x*x.');
% find the second derivative of f - if this is always positive the function
% is strictly convex
ddf = jacobian(jacobian(f,x),x);
% this is now a 3x3 matrix; the equivalent condition for convexity here
% would be a positive definite mtrix
% fTest = eig(ddf);
% or, use the cholesky factorization which returns a positive definite (symmetric) flag
[~,flag] = chol(ddf);
if flag == 1
    fprintf("The function in question 3 is convex!\n")
else
    fprintf("The function in question 3 is nonconvex.\n")
end

%% problem 2.4
fprintf("Q2: \n")
syms alpha c real
assume(c, 'positive')
% bcause the sum of convex functions is itself convex if we can demonstrate
% than any permissible term is convex then the sum must also be convex
g(x) = c*exp(alpha*x1);
ddg = jacobian(jacobian(g,x1),x1);
ddg
fprintf("The function in question 4 is convex; alpha is real so alpha^2 must be positive and c is real. \n Further, e raised to the power of any real number is convex.\n")

%% problem 2.5
fprintf("Q3: \n")
assume([x1,x2],'positive')
h = x1^2 - 4*x1*x2 + 4*x2^2 - log(x1*x2);
ddh = jacobian(jacobian(h,[x1,x2]),[x1,x2]);
[~,flag] = chol(ddh);
if flag == 1
    fprintf("The function in question 5 is convex!\n")
else
    fprintf("The function in question 5 is nonconvex.\n")
end

%% problem 2.6
fprintf("Q4: \n")
assume([x1,x2],'clear')
f = 5*exp(x1) - exp(x2);
ddf = jacobian(jacobian(f,[x1,x2]),[x1,x2]);
eig(ddf)
fprintf("The matrix has one positive and one negative eigenvalue - it is indefinite.\n")

%% problem 2.11
fprintf("Q5: \n")
z = -2:.1:2;
fz = [];
for i = z
    if i < 0
        fz(end+1) = i^2;
    else
        fz(end+1) = i^2 + .75;
    end
end
figure(10101);
plot(z,fz,'LineWidth',4);
title("Plot of Question 11");
fprintf("The function in question 2.11 is quasiconcave do to the discontinuity around the origin.\n")
%% problem 2.17
assume([x1,x2],'positive')
fprintf("Q6: \n")
g = x1*x2;
