function [x,y] = curvatureMap(q,s)
%given state vector q and vector of normalized distances along member (s in
%[0,1]), return x and y coordinates of points at s norm'd dist along member

q0 = q(1); tht = q(2);
[~,L,~,~,~,~,~,~,~,~,~] = getVars();

t2 = q0.*s;
t3 = 1.0./q0;
t4 = t2+tht;
t5 = cos(t4);
x = L.*t3.*(t5-cos(tht));%+D.*d.*t5;
if nargout > 1
    t6 = sin(t4);
    y = L.*t3.*(t6-sin(tht));%+D.*d.*t6;
end
