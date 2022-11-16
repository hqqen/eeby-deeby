clear all; close all;
%% acrobot dynamics simulator
%alex beyer 2022
%this file is contains logic to run the integrator and plot results, the
%control law itself is contained within fdbkCtrl.m, the dynamics are
%generated via a matlabFunction call from acrobot.m (with a custom header
%added later). getSimParams is a helper function which passes paramters
%such as link mass, link length and gravity to the simulation to avoid
%hardcoding them. all other functions miscelaneous plotting helpers

q = [-.1;-.1];%[(-1^randi([1,2]))*.1*rand(1);(-1^randi([1,2]))*.1*rand(1)];
dq = [0;0];
tspan = 0:.01:10;

x = [q;dq];
[t,y] = ode45(@fdbkCtrl,tspan,x);

%figure(); 
subplot(2,2,1); sgtitle("\tau = 0 Plots")
hold on;
plot(0,pi/2,'bx'); 
plot(y(:,2),y(:,1)); legend("Goal","State");
xlabel("\theta_2"); ylabel("\theta_1")
title("State Space Orbit (joint angles)")

%figure(); hold on;
% plot(0,pi/2,'bx');
subplot(2,2,2)
plot(y(:,4),y(:,3));
xlabel("d\theta_2/dt"); ylabel("d\theta_1/dt")
title("State Space Orbit (joint velocities)")

E = []; L = []; T = []; U = [];
for i = 1:max(size(y))
    [h,l,t2,u] = systemEnergy(y(i,:));
    E = [E h]; L = [L l]; T = [T t2]; U = [U u];
end

%figure(); 
subplot(2,2,3)
hold on;
plot(t,abs(E)); plot(t,L);
title("Magnitude of System Energy");
xlabel("Time (s)"); ylabel("Energy (J)")
legend("Hamiltonian","Lagrangian","Location","best")

p1 = []; p2 = [];
for i = 1:max(size(y))
    [p1t,p2t] = getJointPosition(y(i,:));
    p1 = [p1 p1t]; p2 = [p2 p2t];
end

%figure(); 
subplot(2,2,4);
hold on; 
title("End Effector Position Plot");
plot(0,2,'rx'); 
xlim([-2.2 2.2]); ylim([-2.2 2.2]);
plot(p2(1,:),p2(2,:),'b-'); xlabel("x"); ylabel("y");
legend("Goal","Actual");
figure(); hold on; title("Joint 1 Position Plot");
xlim([-2.2 2.2]); ylim([-2.2 2.2]);
plot(p1(1,:),p1(2,:),'r-'); xlabel("x"); ylabel("y");

figure(); hold on;
plot(t,abs(E));
plot(t,abs(T)); plot(t,abs(U));
title("Magnitude of System Energy");
xlabel("Time (s)"); ylabel("Energy (J)")
legend("|Total Energy|", "|Kinetic Energy|","|Poential Energy|")

% figure(); hold on; title("Control Signal"); ylim([-1.1 1.1]);
% plot(t,y(:,5))

% figure();
% for i = 1:max(size(y))
%     plot([0 p1(1,i)],[0 p1(2,i)],'b-'); hold on;
%     plot([p1(1,i) p2(1,i)],[p1(2,i) p2(2,i)],'r-');
%     plot(p1(1,i),p1(2,i),'k*');
%     plot(p2(1,i),p2(2,i),'k*');
%     xlim([-2 2]); ylim([-2 2]);
%     drawnow; hold off;
% end