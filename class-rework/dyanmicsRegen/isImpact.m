function [h, isTerm, direction] = isImpact(t,x)

    q0 = x(1); tht0 = x(2);
    q1 = x(3); tht1 = x(4);
    dq0 = x(5); dtht0 = x(6);
    dq1 = x(7); dtht1 = x(8);
    
    isTerm = 1;
    direction = 0;

    h = swingFootHeight(q0,q1,tht0,tht1);

end