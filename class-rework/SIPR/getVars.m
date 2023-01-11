function [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints] = getVars()
%helper function; returns state vars for computations elsewhere
%D - arm width
%L - arm length
%g - |gravity|
%m - linkage mass (by section?)
%phi - g field angle from negative vertical, ccw

D = .1;
L = 1;
g = 9.81;
m = 2;
phi = 0;
dt = .01;
beta = .1;
ke = 3;
kd = 1.25;
kp = 1.1;
k = .5;
nPoints = 15;

end