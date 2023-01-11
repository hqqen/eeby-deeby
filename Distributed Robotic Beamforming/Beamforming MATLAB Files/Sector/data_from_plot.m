clc
openfig('pattern_1.fig');
h = findobj(gca,'Type','line');
z=get(h,'Zdata');
mean(z)