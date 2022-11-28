f = @(t,x) -x^(4/3);
x0 = -10:1:10;
tspan = 0:.01:3;

figure(); hold on;
for i = x0
    [t, y] = ode45(f,tspan,i);
    plot(t,y); 
end