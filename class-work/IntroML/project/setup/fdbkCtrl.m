function [dx] = fdbkCtrl(t,x)
    %get inputs
    q = x(1:2); dq = x(3:4);
    [M, C, G] = acrobotDyn(x);
    [m1, m2, l1, l2, I1, I2, g] = getSimParams();

    B = [0;1];%control mapping

    %enforce position and speed limits
    dq(1) = min(max(dq(1),-4*pi),4*pi);
    dq(2) = min(max(dq(2),-4*pi),4*pi);
    
    %def control law
    tau = 0;%randi([0 1]);

    %get state update
    dx = [...
         dq;
         M\(B*tau - C*dq - G.')]; %gravity vector malformed!
end