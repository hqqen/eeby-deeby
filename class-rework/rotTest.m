%make sim objects
close all;
simWorld = World;
simWorld.groundPointsX = -2:12; simWorld.groundPointsY = -0.1*simWorld.groundPointsX;
simWorld.makeWorld();

biped = Robot;
% biped.q = [pi/96; 0; pi/24; -pi]; biped.dq = [1; 1.6; 1; -1.6];
biped.world = simWorld;
biped.B = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
q0 = [pi/12; 0; pi/12; 3.5*pi/4]; dq0 = [0;-pi/4;0;-pi/4];%[1; -1.6; 1; -1.6];

biped.drawSelf([q0;dq0],1); title("Pre-SS Rotation")

q0Star = q0;
q0(1) = q0Star(3);
q0(2) = q0Star(4) + pi;
q0(3) = q0Star(1);
q0(4) = q0Star(2) - pi;


biped.drawSelf([q0;dq0],2); title("Post-SS Rotation")