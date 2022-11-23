dtht = 1e-5;
tht = [0:dtht:2*pi];

delta = 1e-15;
y = mod(tan(tht),delta);
yNZ = find(y);

scatter3(cos(tht(yNZ)),sin(tht(yNZ)),y(yNZ),.15); hold on;
plot3(cos(tht),sin(tht),tht*0,'r--','LineWidth',.1);

title("Rational Solutions of tan(\theta) out to \Delta = " + delta, "d\theta = " + dtht)