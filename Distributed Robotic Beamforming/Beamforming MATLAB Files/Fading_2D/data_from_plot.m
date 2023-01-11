clc
openfig('position.fig');
h = findobj(gcf,'Type','scatter');
xm=get(h,'Xdata');
ym=get(h,'Ydata');
xm = xm';
ym = ym';
