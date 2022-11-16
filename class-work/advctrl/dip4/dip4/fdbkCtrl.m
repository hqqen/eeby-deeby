function dx = fdbkCtrl(t,x)
    %uses manipulator form dynamics to simulate the double pendulum for
    %arbirary two-valued control input tau

    Kp = 1;
    Ki = 0;
    Kd = 6;
    
    %get manipulator form dynamics
    q = x(1:2); dq = x(3:4);
    [M,C,G,T,V] = getDynamics(q,dq);
    [I1,I2,m1,m2,a1,a2,k1,k2] = getVars();

    %MAKE CONTROL LAW HERE (only takes one input - tau1)
    tau = [Kp*q(1) + Kd*dq(1);
            0];
    
    %find the second derivative of the states via the manipulator equation
    %and concatenate them with first joint state time derivs to form dx
    ddq = M\(tau - C*dq - G);
    dx = [dq;ddq];
    
end

