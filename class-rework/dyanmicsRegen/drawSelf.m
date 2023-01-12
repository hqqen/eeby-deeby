function fig = drawSelf(y,f,t)
    %given robot state y, figure axes f (opt.) and current sim time
    %t (opt.)
    %return a picture of the robotcs configuration at time t on
    %axes f
    
    if ~exist('f','var')
        %if not given axes, make a set of them way outside what
        %people would normally use
        f = 1000000;
    end
    
    if ~exist('t','var')
        %if not given time, zero it
        t = 000000;
    end
    
    %unpack inputs & get system vars
    q0 = y(1); tht0 = y(2);
    q1 = y(3); tht1 = y(4);
    q = y(1:4); dq = y(5:8);
    footX = 0; footY = 0;
    %[dxSw,dySw,xSw,ySw] = Robot.getSwingFootPos(q); %verify the robot knows where its foot is
    %[D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();
    D = .1; L = 1; nPoints = 10;

    %find hip and point-along-leg positions
    s = 0:1/nPoints:1; d = 0; %only care about centerline
    swingLegX = footX + D.*d.*cos(q0 + tht0) + (L.*(cos(q0 + tht0) - cos(tht0)))./q0 + (L.*(cos(q1 + tht1) - cos(q1 + tht1 - q1.*s)))./q1 + D.*d.*cos(q1 + tht1 - q1.*s);
    swingLegY = footY + D.*d.*sin(q1 + tht1 - q1.*s) + D.*d.*sin(q0 + tht0) + (L.*(sin(q0 + tht0) - sin(tht0)))./q0 + (L.*(sin(q1 + tht1) - sin(q1 + tht1 - q1.*s)))./q1;
    
    stanceLegX = footX + D.*d.*cos(tht0 + q0.*s) + (L.*(cos(tht0 + q0.*s) - cos(tht0)))./q0;
    stanceLegY = footY + D.*d.*sin(tht0 + q0.*s) + (L.*(sin(tht0 + q0.*s) - sin(tht0)))./q0;
    
    hipX = footX + D.*d.*cos(q0 + tht0) + (L.*(cos(q0 + tht0) - cos(tht0)))./q0;
    hipY = footY + D.*d.*sin(q0 + tht0) + (L.*(sin(q0 + tht0) - sin(tht0)))./q0;
    
    %plot the biped
    figure(f); hold off; %we're animating on this, forcibly wipe it while drawing the biped
    plot([stanceLegX hipX],[stanceLegY hipY],'b-'); hold on;
    plot([-10 10], [0 0], 'k-');
    ylim([hipY - 1.25, hipY + 1.25]); xlim([hipX - 1.25, hipX + 1.25]); %center plot on hip, makin sure it has a window big eough to stop limbs from leaving
    plot([hipX swingLegX], [hipY swingLegY],'r-');
    %plot(xSw,ySw,'k*'); %plot separately calc'd foot position to check for issues
    title("Robot State at " + t + " sec");
    %plot(Robot.world.groundPointsX,Robot.world.groundPointsY); %plot ground
    fig = gcf;
end