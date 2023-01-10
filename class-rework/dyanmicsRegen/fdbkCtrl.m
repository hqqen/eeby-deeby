function [dx] = fdbkCtrl(t,x)

    q0 = x(1); tht0 = x(2);
    q1 = x(3); tht1 = x(4);
    q = [q0; tht0; q1; tht1];
    dq0 = x(5); dtht0 = x(6);
    dq1 = x(7); dtht1 = x(8);
    dq = [dq0; dtht0; dq1; dtht1];

    [M,C,G] = genDynamics(dq0,dq1,dtht0,dtht1,q0,q1,tht0,tht1);
    [T,P] = getEnergy(dq0,dq1,dtht0,dtht1,q0,q1,tht0,tht1);
    %for group action add the slope angle in rads to each base angle
    Gphi = getGroupAction(q0,q1,tht0,tht1);

    k = 1.13; b = .1;
    K = [k 0 0 0;...
         0 0 0 0;
         0 0 k 0;
         0 0 0 0];
    B = [b 0 0 0;
         0 0 0 0;
         0 0 b 0;
         0 0 0 0];

    %find control input tau:
    E = T + P + .5*q.'*K*q;
    Eref = P + 1;
    kp = .1;
    tau = eye(4)\(G.' - Gphi.' - kp*eye(4)\((E-Eref)*dq));
    tau(2) = 0;

    %should have inputs on 1st and 3rd indices but only actuating base
    %angle for now to test rigid gait generation
    dx = M\(tau  - C*dq - B*dq - K*q - G.');

    %zero out curvature RoCs
    dx = [dq; dx];
    dx(1) = 0; dx(3) = 0; dx(5) = 0; dx(7) = 0;
    
end

