syms s l real
L = 3; tht0 = pi/24; q0 = pi/12

alpha0 = tht0 + s*q0;

xsd0 = -L*int(sin(alpha0),s,0,s);
ysd0 = L*int(cos(alpha0),s,0,s);

xbar = L*int(-sin(alpha0),s,0,1);
ybar = L*int(cos(alpha0),s,0,1);

%find model CoM (draw a line connecting the two ends of the link and place
%CoM on the exact center)
xsd0i = subs(xsd0,s,0); ysd0i = subs(ysd0,s,0);
xsd0f = subs(xsd0,s,1); ysd0f = subs(ysd0,s,1);
projX = .5*(xsd0i + xsd0f); 
projY = .5*(ysd0i + ysd0f);

%plot
s = 0:.1:1;
memberX = (L.*(cos(tht0 + q0.*s) - cos(tht0)))./q0;
memberY = (L.*(sin(tht0 + q0.*s) - sin(tht0)))./q0;

plot(memberX, memberY, 'k-', 'LineWidth',2); hold on;
plot(xbar, ybar, 'b*')
plot([xsd0i, xsd0f],[ysd0i, ysd0f],'r-')
plot(projX, projY, 'r*')