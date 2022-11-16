function [I1,I2,m1,m2,a1,a2,k1,k2] = getVars()
    %helper function to return system parameters to double inverted
    %pendulum integrator

    m1 = 10;
    m2 = 20;
    a1 = 1;
    a2 = 1;
    k1 = 5;
    k2 = 30;
    I1 = m1*a1^2; %this will always be correct
    I2 = m2*a2^2; %this gets less correct the more the 2nd joint bends

end