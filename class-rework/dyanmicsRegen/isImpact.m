function [h, isTerm, direction] = isImpact(t,x)

    q0 = x(1); tht0 = x(2);
    q1 = x(3); tht1 = x(4);
    
    isTerm = 1;
    direction = 1;

    h = swingFootHeight(q0,q1,tht0,tht1);

end