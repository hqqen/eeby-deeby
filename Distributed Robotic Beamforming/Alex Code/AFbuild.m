close all
tht = 0:pi/100:2*pi;

for i = 2:max(size(tht))

    y = cos(tht) + sin(tht);

end

plot(tht,y)