sys = ss(A2,B2,C2,D2);
Q = eye(4);
fprintf("LQR Gains Found to BE:\n");
[K, S, P] = lqr(sys,Q,1);
disp(lqr(sys,Q,1))
N = 10;

stab_sys = ss(A2-B2*K,B2*N,C2,D2);

%load in the linearized system matrices generated by DIPdyn.m
load("getLinearizedSystem.mat",'-mat')

%this will give 4 transfer functions per input, we only care about the
%first 2! (the second two are joint vvelocities which we dont want)
stab_pend = tf(ss(A2-B2*K,B2,C2,D2)); %stab_stab_pend(i,j) is the jth variables effect on the ith ouput
size(stab_pend(1,1))
fb = bandwidth(stab_pend(1,1));
disp("Bandwidth = ")
disp(fb)

%figure();
step(stab_pend(1,1))
disp(stepinfo(stab_pend(1,1)))

figure();
sgtitle("Stabilized Root Loci");
subplot(2,1,1)
rlocus(stab_pend(1,1));
title("\tau_1 -> \theta_1");
subplot(2,1,2);
rlocus(stab_pend(2,1));
title("\tau_1 -> \theta_2");

figure();
bode(stab_pend(1,1));
title("Stabilized Bode Plot (\tau_1 -> \theta_1)");
figure();
bode(stab_pend(2,1));
title("Stabilized Bode Plot (\tau_1 -> \theta_2)");

figure();
sgtitle("Stabilized Nyquist Plots")
subplot(2,1,1)
nyquist(stab_pend(1,1))
title("\tau_1 -> \theta_1");
subplot(2,1,2)
nyquist(stab_pend(2,1))
title("\tau_1 -> \theta_2");


autoArrangeFigures(2,2,1)


