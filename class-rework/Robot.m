classdef Robot
    %UNTITLED Summary of this class goes here
    %   Detailedexplanation goes here

    properties
        %initialize an at rest robot
        %         q = [0;0;0;0]; %[qSt, thtSt, qSw, thtSw]
        %         dq = [0;0;0;0]; %[dqSt, dthtSt, dqSw, dthtSw]
        p = [0;0]; %want to initialize swing foot for any arbitrary state @ origin
        %x = [q;dq]; %unneeded; state gets passed out by methods which will build the outgoing state OTF
        dynFlag = 0;
        isFix = 0; %flag to check if the inputs to the robot have been columnized
        mode = 1; %what controller to use
        world = ""; %what worldspace the robot is operating in
        B %input mapping for running feedback control
        R = eye(4);
        tol = 1e-3;
    end

    methods

        %         function [M,C,G,E,PE,KE] = doStep(Robot,q,dq)
        %             %find xdot (dx) for current step by solving the modified E-L
        %             %equation, also functions as a dynamics listener if more than
        %             %one output is specified
        %             %returns dx
        %             %want to take in the q and dq (given as x) from integrator, but
        %             %use the stance foot posn of the robot due to how ode45 works
        %             if Robot.dynFlag == 0
        %                 if nargout == 2
        %                     fprintf("Listening to Controller Input...")
        %                     Robot.dynFlag = 1;
        %                 elseif nargout > 2
        %                     fprintf("Listening to Dynamics...")
        %                     Robot.dynFlag = 1;
        %                 else
        %                     fprintf("Only Reporting dx...")
        %                     Robot.dynFlag = 1;
        %                 end
        %             end
        %             if Robot.isFix == 0
        %                 Robot.isFix = 1;
        %
        %             end
        %             [M,C,G,E,PE,KE] = Robot.getDyns(q,dq,Robot.p);
        %         end

        function h = getHeight(Robot,q)
            [~,~,xSw,ySw] = Robot.getSwingFootPos(q); %use integrator made q for this
            h = Robot.world.getHeight(xSw,ySw);
        end

        function [h, isterminal, direction] = isImpact(Robot,t,y)
            isterminal = 1;
            direction = -1;
            %add condition here to check for swing foot being right of
            %stance foot (forcibly stop rocking)
            [~,~,xSw,ySw] = Robot.getSwingFootPos(y);
            if xSw > Robot.p(1)
                h = ySw;
            else
                h = inf;
            end
            %h = (Robot.getHeight(y(1:4)) < 0); %ode45 can't interpolate a function call like this
        end

        function y = doImpact(Robot,q,dq)
            %run impact dynamics on current robot state
            n = max(size(q));
            [Me,E2,dGamma] = Robot.getImpactDyns(q,dq,Robot.p);
            [~,~,xSw,ySw] = Robot.getSwingFootPos(q);
            %             R = [0 0 1 -1;...
            %                  0 0 0 1;
            %                  1 1 0 0;
            %                  0 1 0 0];
            %             R = [0 0 -1 0;...
            %                  0 0 0 -1;
            %                  -1 0 0 0;
            %                  0 -1 0 0];

            dqPlus = inv(-E2*inv(Me)*E2.');
            dqPlus = (Me\(E2.'*dqPlus));
            dqPlus = dqPlus*E2*[eye(n);dGamma];
            dqPlus = dqPlus + [eye(n);dGamma];
            dqPlus = [Robot.R, zeros(n,2)]*dqPlus*dq.';

            Robot.p = [xSw;ySw];
            xPlus = [Robot.R*q.'; dqPlus]; %put foot coords back in

            %have the new state, need to rearrange:
            yStar = xPlus;
            y = [...
                yStar(3); %new stance leg curvature is old swing leg curvature
                yStar(4) + pi; %new stance leg base angle is old swing leg base angle w/offset
                yStar(1); %new swing leg curvature is old stance leg curvature
                yStar(2) - pi; %new swing leg base angle is old stance leg base angle w/offset
                yStar(7); %new stance leg curvature RoC is old swing leg curvature RoC
                yStar(8); %new stance leg base angle RoC is old swing leg base angle RoC w/o an offset
                yStar(5); %new swing leg curvature RoC is old stance leg curavture RoC
                yStar(6)]; %new swing leg base angle RoC is old stance leg base angle RoC w/o an offset
            %             y = xPlus;
        end

        function dx = doCtrl(Robot,t,x)
            q = x(1:4); dq = x(5:8);

            %singularities occur for straightened legs, set a tolerance and
            %round things away from 0
            tol = Robot.tol;
            if abs(q(1)) < tol
                if dq(1) < 0
                    q(1) = -tol;
                else
                    q(1) = tol;
                end
            end
            if abs(q(3)) < tol
                if dq(3) < 0
                    q(3) = -tol;
                else
                    q(3) = tol;
                end
            end


            %picks frm a library of conrollers to get tau
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();
            [M,C,G,~,~,~] = Robot.getDynsSpring(q,dq,Robot.p);
            tau = Robot.getTau(M,C,G,q,dq);

            %member stiffness and damping apply to both leg curvatures
            ddq = M\(tau - C*dq + [beta*dq(1); 0; beta*dq(3); 0] + [k*q(1); 0; k*q(3); 0] + G); %shouldnt G be negative here? ----%#%@%$^
            dx = [dq;ddq];
        end

        function tau = getTau(Robot,M,C,G,q,dq)
            %picks frm a library of conrollers to get tau
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();
            qDes = pi/24; %put in getvars
            %build control input
            %build slope correcting group action
            %only applies to base angles (doing this to curvatures would
            %result in permenently curved members)
            Gphi = Robot.doGroupAction(q,dq);
            tauGrav = pinv(Robot.B)*(G-Gphi); %slope correction
            %build straighteneing action
            tauStraighten = [-kd*dq(1) - kp*(q(1)-qDes);0; -kd*dq(3) - kp*(q(1)+qDes);0]; %note that q3 is given the oppsotie qDes from q1
            %build rigid gait input
            %we've corrected for slope and straightened our legs, use
            %motion at the hip to generate the gait itself
            [~, lffh, lglfh, ~] = Robot.getLieDeriv(q,dq,M,C,G);
            tauHip = lglfh\lffh;%lglfh\(Robot.getPsi(x(1:4),x(5:8)) - lffh);
            tau = tauGrav + tauStraighten + tauHip;
        end

        function [lfh, lffh, lglfh, psi] = getLieDeriv(Robot,q,dq,M,C,G)
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            %TODO: convert into function call (not done here b/c objective
            %function is simple enough to turn both derivatives into
            %constant matrices)
            [h,~,dhdq,d2hdq2] = Robot.getObjFunc(q,dq);

            lfh = [dhdq]*[dq.'; (M\(-C*dq - G)).'].'; %don't need the zeros as in wvelt ??
            lffh = [d2hdq2 dhdq]*[dq; M\(-C*dq - G)];
            lglfh = [dhdq]*[M\Robot.B];

            psi = Robot.getPsi(h,lfh);%dh <--> lfh

        end

        function [h,dh,dhdq,d2hdq2] = getObjFunc(Robot,q,dq)
            %calculates objective function h for tracking
            %additionally returns first time derivative dhdt and (1st & 2nd)
            %state derivatives dhdq and d2hdq2
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            %find objective function at current state
            h = [q(1) - tht0des;...
                q(1) - q(3)];

            %find first and second derivs w/r.t. state vars
            %constant matrices for now, but will get more complex as
            %objective function does
            dhdq = [1, 0, 0, 0;...
                1, 0, -1, 0];
            d2hdq2 = [0, 0, 0, 0;...
                0, 0, 0, 0];

            %find first time deriv
            dh = [dq(1);...
                 dq(1) + dq(3)];
        end

        function [psi] = getPsi(Robot,x,dx)
            %get psi function for hip controller
            %taken as in grizzle, '01 but is selected from famly of valid
            %ctrlrs
            x1 = x(1); x2 = x(2); dx1 = dx(1); dx2 = dx(2);
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            psi = (1/epsilon^2) * [...
                x1 + (-alpha + .5)*sign(epsilon*dx1)*abs(epsilon*dx1)^(2-alpha);
                x2 + (-alpha + .5)*sign(epsilon*dx2)*abs(epsilon*dx2)^(2-alpha)];

        end

        function fig = animate(Robot,y,t,f,fixed)
            %given a series of states y evolving w/r.t. time vector t
            %(opt.)
            %(max(size(y)) == size(t)) plot the physical state of teh robot
            %on figure axes f (opt.) and animate them
            %also takes in a toggle named fixed (opt.), 0 means camera
            %centered at hip, 1 means centered at (0,0)

            if ~exist('f','var')
                %if not given axes, make a set of them way outside what
                %people would normally use
                f = 1000000;
            end

            if ~exist('t','var')
                %if not given time, zero it
                t = 000000;
            end

            if ~exist('fixed','var')
                %if not given fix state, assume centering at hip
                fixed = 000000;
            end

            %animate by repeatedly calling the drawSelf function, fixing camera if
            %necessary
            if fixed == 0
                for i = 1:max(size(y))
                    fig = Robot.drawSelf(y(i,:),f,t(i)); drawnow;
                end
            else
                for i = 1:max(size(y))
                    fig = Robot.drawSelfFixedCam(y(i,:),f,t(i)); drawnow;
                end
            end
        end

        function fig = drawSelf(Robot,y,f,t)
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
            footX = Robot.p(1); footY = Robot.p(2);
            [dxSw,dySw,xSw,ySw] = Robot.getSwingFootPos(q); %verify the robot knows where its foot is
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

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
            ylim([hipY - 1, hipY + 1]); xlim([hipX - 1, hipX + 1]); %center plot on hip, makin sure it has a window big eough to stop limbs from leaving
            plot([hipX swingLegX], [hipY swingLegY],'r-');
            plot(xSw,ySw,'k*'); %plot separately calc'd foot position to check for issues
            title("Robot State at " + t + " sec");
            plot(Robot.world.groundPointsX,Robot.world.groundPointsY); %plot ground
            fig = gcf;
        end

        function fig = drawSelfFixedCam(Robot,y,f,t)
            %given robot state y, figure axes f (opt.) and current sim time
            %t (opt.)
            %return a picture of the robotcs configuration at time t on
            %axes f
            %this version of the function fixes the robot w/in the window
            %(in particular, the stance foot at (0,0)) to stop weid perspective effects from changing how it's moving

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
            footX = 0; footY = 0; %lock foot pos'n
            [dxSw,dySw,xSw,ySw] = Robot.getSwingFootPos(q); %verify the robot knows where its foot is
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            %find hip and point-along-leg positions
            s = 0:1/nPoints:1; d = 0; %only care about centerline
            swingLegX = footX + D.*d.*cos(tht1 + q1.*s) + (L.*(cos(tht1 + q1.*s) - cos(tht1)))./q1 + D.*d.*cos(q0 + tht0) + (L.*(cos(q0 + tht0) - cos(tht0)))./q0;
            swingLegY = footY + D.*d.*sin(tht1 + q1.*s) + (L.*(sin(tht1 + q1.*s) - sin(tht1)))./q1 + D.*d.*sin(q0 + tht0) + (L.*(sin(q0 + tht0) - sin(tht0)))./q0;

            stanceLegX = footX + D.*d.*cos(tht0 + q0*s) + (L.*(cos(tht0 + q0.*s) - cos(tht0)))./q0;
            stanceLegY = footY + D.*d.*sin(tht0 + q0*s) + (L.*(sin(tht0 + q0.*s) - sin(tht0)))./q0;

            hipX = footX + D.*d.*cos(q0 + tht0) + (L.*(cos(q0 + tht0) - cos(tht0)))./q0;
            hipY = footY + D.*d.*sin(q0 + tht0) + (L.*(sin(q0 + tht0) - sin(tht0)))./q0;

            %plot the biped
            figure(f); hold off; %we're animating on this, forcibly wipe it while drawing the biped
            plot([stanceLegX hipX],[stanceLegY hipY],'b-'); hold on;
            ylim([-1 1]); xlim([-1 1]); %center plot on hip, makin sure it has a window big eough to stop limbs from leaving
            plot([hipX swingLegX], [hipY swingLegY],'r-');
            plot(xSw,ySw,'k*'); %plot separately calc'd foot position to check for issues
            title("Robot State at " + t + " sec");
            plot(Robot.world.groundPointsX,Robot.world.groundPointsY); %plot ground
            fig = gcf;
        end

        function [Gphi] = doGroupAction(Robot,q,dq)
            [~,~,xSw,~] = getSwingFootPos(Robot,q);
            phi = Robot.world.getSlope(xSw);
            qPhi = q;
            qPhi(1) = qPhi(1) + phi; qPhi(3) = qPhi(3) + phi; %should onyl change base angles, not curavtures
            [~,~,Gphi,~,~,~] = Robot.getDynsSpring(qPhi,dq,Robot.p); %w/ or w/o spring in PE? --tested and doesnt seem to matter
        end

        function [dxSw,dySw,xSw,ySw] = getSwingFootPos(Robot,q)
            %helper function to return swing foot velocity and position for impact detection

            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();
            q0 = q(1); tht0 = q(2); q1 = q(3); tht1 = q(4);
            footX = Robot.p(1); footY = Robot.p(2); d = 0;

            t2 = cos(tht0);
            t3 = cos(tht1);
            t4 = sin(tht0);
            t5 = sin(tht1);
            t6 = q0+tht0;
            t7 = q1+tht1;
            t12 = 1.0./q0;
            t14 = 1.0./q1;
            t8 = cos(t6);
            t9 = cos(t7);
            t10 = sin(t6);
            t11 = sin(t7);
            t13 = t12.^2;
            t15 = t14.^2;
            t16 = -t2;
            t17 = -t3;
            t18 = -t4;
            t19 = -t5;
            t20 = D.*d.*t3;
            t21 = D.*d.*t5;
            t22 = D.*d.*t8;
            t23 = D.*d.*t10;
            t25 = t8+t16;
            t26 = t9+t17;
            t27 = t10+t18;
            t28 = t11+t19;
            t29 = -L.*t12.*(t2-t8);
            t30 = -L.*t14.*(t3-t9);
            t31 = -L.*t12.*(t4-t10);
            t32 = -L.*t14.*(t5-t11);
            t24 = -t23;
            dxSw = [t24-L.*t10.*t12+L.*t13.*(t2-t8),t24+L.*t12.*(t4-t10),-L.*t11.*t14+L.*t15.*(t3-t9),D.*d.*t19+L.*t14.*(t5-t11)];
            if nargout > 1
                dySw = [t22+L.*t8.*t12+L.*t13.*(t4-t10),t22+t29,L.*t9.*t14+L.*t15.*(t5-t11),t20+t30];
            end
            if nargout > 2
                xSw = footX+t20+t22+t29+t30;
            end
            if nargout > 3
                ySw = footY+t21+t23+t31+t32;
            end
        end

        function [Me,E2,dGamma] = getImpactDyns(Robot,q,dq,p)
            %helper function to numerically calculate extended inertia matrix Me for
            %impact mapping

            q0 = q(1); tht0 = q(2); q1 = q(3); tht1 = q(4);
            dq0 = dq(1); dtht0 = dq(2); dq1 = dq(3); dtht1 = dq(4);
            footX = p(1); footY = p(2); d = 0;

            % fprintf("q:\n"); disp(q);
            % fprintf("\ndq:\n"); disp(dq);
            % fprintf("\nqo:\n"); disp(q0);
            % fprintf("\ntht:\n"); disp(tht);
            % fprintf("\ndq0:\n"); disp(dq0);
            % fprintf("\ndtht:\n"); disp(dtht);

            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars(); d = 0;

            t2 = cos(q0);
            t3 = cos(q1);
            t4 = cos(tht0);
            t5 = cos(tht1);
            t6 = sin(q0);
            t7 = sin(q1);
            t8 = sin(tht0);
            t9 = sin(tht1);
            t10 = q0+tht0;
            t11 = q1+tht1;
            t12 = D.^2;
            t13 = L.^2;
            t14 = m.*3.0;
            t15 = q0.^2;
            t16 = q0.^3;
            t17 = q1.^2;
            t19 = q1.^3;
            t29 = -q1;
            t30 = 1.0./q0;
            t32 = 1.0./q1;
            t38 = -tht0;
            t39 = -tht1;
            t18 = t15.^2;
            t20 = t17.^2;
            t21 = t4.*2.0;
            t22 = t5.*2.0;
            t23 = t8.*2.0;
            t24 = t9.*2.0;
            t25 = cos(t10);
            t26 = cos(t11);
            t27 = sin(t10);
            t28 = sin(t11);
            t31 = 1.0./t15;
            t33 = 1.0./t16;
            t34 = 1.0./t17;
            t36 = 1.0./t19;
            t40 = -t4;
            t42 = -t5;
            t44 = -t8;
            t50 = t13.*8.0;
            t51 = t13.*2.4e+1;
            t58 = t39+tht0;
            t59 = q0.*t4.*-2.0;
            t64 = t10+t39;
            t65 = t11+t38;
            t87 = t13.*t15.*1.2e+1;
            t88 = t13.*t17.*3.6e+1;
            t89 = t2.*t13.*-8.0;
            t93 = q0.*t6.*t13.*-8.0;
            t94 = q1.*t7.*t13.*-2.4e+1;
            t97 = t10-t11;
            t35 = 1.0./t18;
            t37 = 1.0./t20;
            t41 = -t21;
            t43 = -t22;
            t45 = -t23;
            t46 = q0.*t21;
            t47 = q0.*t23;
            t48 = t25.*2.0;
            t49 = t26.*2.0;
            t52 = t27.*2.0;
            t53 = t28.*2.0;
            t54 = q0.*t25;
            t55 = q1.*t26;
            t56 = q0.*t27;
            t57 = q1.*t28;
            t60 = -t51;
            t61 = cos(t58);
            t62 = -t28;
            t66 = sin(t58);
            t69 = D.*d.*t25;
            t70 = D.*d.*t27;
            t71 = t12.*t18;
            t72 = t12.*t20;
            t73 = cos(t64);
            t74 = cos(t65);
            t76 = sin(t64);
            t77 = sin(t65);
            t78 = t15.*t25;
            t79 = t17.*t26;
            t80 = t2.*t50;
            t81 = t3.*t51;
            t82 = t15.*t27;
            t83 = t17.*t28;
            t84 = q0.*t6.*t50;
            t85 = q1.*t7.*t51;
            t90 = L.*t25.*t30;
            t91 = t25+t40;
            t92 = t27+t44;
            t98 = cos(t97);
            t99 = sin(t97);
            t100 = t15.*t89;
            t126 = -L.*t30.*(t4-t25);
            t130 = -L.*t31.*(t8-t27);
            t143 = L.*t31.*(t8-t27);
            t63 = -t53;
            t67 = q1.*t49;
            t68 = q1.*t53;
            t75 = -t56;
            t86 = -t70;
            t95 = -t79;
            t96 = t15.*t80;
            t101 = t13.*t61.*1.2e+1;
            t102 = t51.*t61;
            t103 = t13.*t66.*1.2e+1;
            t105 = t13.*t61.*-2.4e+1;
            t106 = q0.*q1.*t12.*t73;
            t107 = q0.*q1.*t12.*t76;
            t108 = t13.*t73.*1.2e+1;
            t109 = t13.*t74.*1.2e+1;
            t110 = t51.*t73;
            t111 = t51.*t74;
            t112 = m.*t51.*t66;
            t113 = t13.*t76.*1.2e+1;
            t114 = t13.*t77.*1.2e+1;
            t115 = m.*t51.*t76;
            t116 = m.*t51.*t77;
            t119 = q1.*t51.*t77;
            t122 = m.*t13.*t76.*-2.4e+1;
            t123 = q0.*t12.*t29.*t73;
            t127 = t26+t42+t57;
            t128 = q1.*t12.*t15.*t76;
            t129 = t9+t55+t62;
            t131 = q1.*t13.*t74.*-1.2e+1;
            t132 = m.*q1.*t12.*t15.*t73;
            t134 = m.*q1.*t13.*t74.*-2.4e+1;
            t135 = q0.*q1.*t12.*t98;
            t136 = q0.*q1.*t12.*t99;
            t137 = t13.*t98.*1.2e+1;
            t138 = t51.*t98;
            t139 = t13.*t17.*t74.*-1.2e+1;
            t140 = t13.*t99.*1.2e+1;
            t141 = t12.*t15.*t29.*t76;
            t144 = m.*t12.*t15.*t29.*t73;
            t146 = t13.*t98.*-2.4e+1;
            t148 = m.*t51.*t99;
            t149 = m.*t13.*t17.*t77.*-1.2e+1;
            t152 = q1.*t51.*t99;
            t156 = q0.*t12.*t29.*t99;
            t157 = q0.*t12.*t17.*t99;
            t158 = q1.*t12.*t15.*t99;
            t159 = q0.*t13.*t99.*-1.2e+1;
            t162 = q0.*q1.*t13.*t98.*-1.2e+1;
            t164 = m.*q1.*t12.*t15.*t98;
            t165 = t54+t59+t92;
            t169 = m.*t12.*t15.*t17.*t99;
            t171 = m.*t13.*t17.*t99.*-1.2e+1;
            t174 = t69+t126;
            t175 = t41+t47+t48+t78;
            t178 = t45+t52+t59+t82;
            t180 = -L.*m.*t31.*(t8-t27+t46-t54);
            t182 = -L.*m.*t31.*(t4-t25+t56-q0.*t8.*2.0);
            t184 = -L.*m.*t33.*(t23+t46-t52-t82);
            t186 = L.*m.*t33.*(t23+t46-t52-t82);
            t188 = t69+t90+t143;
            t189 = -L.*m.*t36.*(t22-t49-t57.*2.0+t79);
            t190 = L.*m.*t36.*(t22-t49-t57.*2.0+t79);
            t191 = t60+t72+t81+t88+t94;
            t193 = t50+t71+t87+t89+t93+t100;
            t104 = -t101;
            t117 = q0.*t113;
            t118 = q1.*t114;
            t120 = -t113;
            t121 = q1.*t109;
            t124 = m.*q0.*t110;
            t125 = m.*q1.*t111;
            t133 = t17.*t109;
            t142 = m.*t17.*t114;
            t145 = -t137;
            t147 = q1.*t137;
            t150 = q0.*t140;
            t151 = q1.*t140;
            t153 = m.*q0.*t138;
            t154 = m.*q1.*t138;
            t160 = t17.*t137;
            t161 = m.*q0.*t146;
            t163 = q0.*q1.*t148;
            t166 = -t157;
            t167 = m.*t17.*t140;
            t170 = t47+t75+t91;
            t172 = L.*m.*t34.*t127;
            t173 = L.*m.*t34.*t129;
            t176 = -t169;
            t179 = t24+t63+t67+t83;
            t181 = t43+t49+t68+t95;
            t183 = L.*m.*t33.*t175;
            t192 = (m.*t37.*t191)./2.4e+1;
            t194 = (m.*t35.*t193)./8.0;
            t155 = q0.*t147;
            t168 = m.*q0.*t160;
            t177 = -t172;
            t185 = L.*m.*t36.*t179;
            t195 = t103+t107+t114+t120+t131+t140+t147+t156;
            t197 = t104+t108+t109+t117+t118+t141+t145+t151+t158+t159+t162;
            t198 = t105+t110+t111+t119+t123+t135+t139+t146+t152+t160+t166;
            t187 = -t185;
            t196 = (m.*t30.*t34.*t195)./1.2e+1;
            t199 = (m.*t31.*t34.*t197)./1.2e+1;
            t201 = (m.*t30.*t36.*t198)./1.2e+1;
            t202 = t112+t116+t122+t124+t134+t144+t148+t149+t154+t161+t163+t164+t168+t171+t176;
            t200 = -t199;
            t203 = (t31.*t36.*t202)./1.2e+1;
            Me = reshape([(m.*t30.^5.*(q0.*t13.*3.6e+1-t6.*t13.*3.6e+1+t13.*t16.*1.2e+1+t12.*1.0./t30.^5-t6.*t13.*t15.*1.8e+1))./9.0,t194,t203,t200,t186,t183,t194,(m.*t33.*(q0.*t51-t6.*t13.*1.2e+1+t12.*t16-q0.*t2.*t13.*1.2e+1))./6.0,t201,t196,t182,t180,t203,t201,(m.*t32.^5.*(t7.*t13.*-1.44e+2+t13.*t19.*4.8e+1+t12.*1.0./t32.^5+q1.*t3.*t13.*1.44e+2))./3.6e+1,t192,t187,t190,t200,t196,t192,(m.*t36.*(q1.*t51-t7.*t13.*2.4e+1+t12.*t19))./1.2e+1,t177,t173,t186,t182,t187,t177,t14,0.0,t183,t180,t190,t173,0.0,t14],[6,6]);
            if nargout > 1
                E2 = [t188,t174,L.*t34.*(t9+t62)+L.*t26.*t32,D.*d.*t5-L.*t32.*(t5-t26),0.0,1.0];
            end
            if nargout > 2
                dGamma = reshape([t86-L.*t27.*t30+L.*t31.*(t4-t25),t188,t86+L.*t30.*(t8-t27),t174,0.0,0.0,0.0,0.0],[2,4]);
            end
        end

        function [M,C,G,E,PE,KE] = getDyns(Robot,q,dq,p)
            %helper function to numerically calculate dynamics matrices M (inertia), C
            %(coriolis) and G (gravity) given state vector q and its derivatives along
            %with hip pos'n vector p

            q0 = q(1); tht0 = q(2); q1 = q(3); tht1 = q(4);
            dq0 = dq(1); dtht0 = dq(2); dq1 = dq(3); dtht1 = dq(4);
            footX = p(1); footY = p(2);


            % fprintf("q:\n"); disp(q);
            % fprintf("\ndq:\n"); disp(dq);
            % fprintf("\nqo:\n"); disp(q0);
            % fprintf("\ntht:\n"); disp(tht);
            % fprintf("\ndq0:\n"); disp(dq0);
            % fprintf("\ndtht:\n"); disp(dtht);

            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            t2 = cos(phi);
            t3 = cos(q0);
            t4 = cos(q1);
            t5 = sin(phi);
            t6 = sin(q0);
            t7 = sin(q1);
            t8 = phi+tht0;
            t9 = phi+tht1;
            t10 = D.^2;
            t11 = L.^2;
            t12 = q1.*2.0;
            t13 = q0.*4.0;
            t14 = q0.^2;
            t15 = q0.^3;
            t16 = q1.^2;
            t18 = q1.^3;
            t19 = q0.^5;
            t21 = q1.^5;
            t31 = -q1;
            t33 = 1.0./q0;
            t35 = 1.0./q1;
            t43 = -tht0;
            t44 = -tht1;
            t17 = t14.^2;
            t20 = t16.^2;
            t22 = t6.*3.0;
            t23 = t7.*3.0;
            t24 = q0.*t3;
            t25 = q1.*t4;
            t26 = cos(t8);
            t27 = cos(t9);
            t28 = q0+t8;
            t29 = q1+t9;
            t30 = sin(t8);
            t32 = -t13;
            t34 = 1.0./t14;
            t36 = 1.0./t15;
            t37 = 1.0./t16;
            t39 = 1.0./t18;
            t40 = 1.0./t19;
            t42 = 1.0./t21;
            t50 = t11.*8.0;
            t51 = t11.*2.4e+1;
            t52 = t6.*t14;
            t53 = t44+tht0;
            t59 = q1+t43+tht1;
            t61 = footY.*g.*m.*t2.*2.0;
            t64 = q0.*t11.*3.6e+1;
            t65 = footX.*g.*m.*t5.*2.0;
            t66 = t10.*t15;
            t68 = t10.*t18;
            t69 = t10.*t19;
            t71 = t10.*t21;
            t76 = (k.*t14)./2.0;
            t77 = (k.*t16)./2.0;
            t80 = t6.*t11.*1.2e+1;
            t82 = t6.*t11.*3.6e+1;
            t86 = t11.*t14.*1.2e+1;
            t87 = t11.*t15.*1.2e+1;
            t88 = t11.*t16.*3.6e+1;
            t89 = t11.*t18.*4.8e+1;
            t90 = t3.*t11.*-8.0;
            t92 = t7.*t11.*-2.4e+1;
            t94 = t7.*t11.*1.44e+2;
            t97 = q0.*t6.*t11.*-8.0;
            t38 = 1.0./t17;
            t41 = 1.0./t20;
            t45 = cos(t28);
            t46 = cos(t29);
            t47 = sin(t28);
            t48 = sin(t29);
            t49 = -t23;
            t54 = -t26;
            t55 = -t27;
            t56 = -t51;
            t57 = cos(t53);
            t58 = q0+t53;
            t60 = sin(t53);
            t62 = q0.*t51;
            t63 = q1.*t51;
            t67 = t10.*t17;
            t70 = t10.*t20;
            t73 = cos(t59);
            t75 = sin(t59);
            t78 = t3.*t50;
            t79 = t4.*t51;
            t81 = t7.*t51;
            t83 = t11.*t24.*1.2e+1;
            t84 = q0.*t6.*t50;
            t91 = -t80;
            t93 = -t82;
            t96 = t11.*t25.*1.44e+2;
            t98 = q1.*t92;
            t99 = -t94;
            t101 = t11.*t52.*1.8e+1;
            t105 = t14.*t90;
            t107 = L.*g.*m.*t30.*t33.*2.0;
            t159 = t22+t24+t32+t52;
            t72 = cos(t58);
            t74 = sin(t58);
            t85 = t7.*t63;
            t95 = -t83;
            t100 = t14.*t78;
            t102 = t31+t58;
            t106 = -t101;
            t108 = t11.*t57.*1.2e+1;
            t109 = t51.*t57;
            t110 = t12+t25+t49;
            t111 = t11.*t60.*1.2e+1;
            t112 = L.*g.*m.*t33.*t45;
            t113 = L.*g.*m.*t35.*t46;
            t114 = L.*g.*m.*t33.*t47;
            t115 = L.*g.*m.*t35.*t48;
            t116 = t45+t54;
            t117 = t46+t55;
            t118 = -t107;
            t120 = t11.*t57.*-2.4e+1;
            t124 = t11.*t73.*1.2e+1;
            t126 = t51.*t73;
            t127 = m.*t51.*t60;
            t129 = t11.*t75.*1.2e+1;
            t131 = m.*t51.*t75;
            t134 = t63.*t75;
            t143 = m.*t63.*t73;
            t145 = q1.*t11.*t73.*-1.2e+1;
            t148 = m.*q1.*t11.*t73.*-2.4e+1;
            t155 = t11.*t16.*t73.*-1.2e+1;
            t165 = m.*t11.*t16.*t75.*-1.2e+1;
            t183 = -L.*g.*m.*t34.*(t26-t45);
            t184 = -L.*g.*m.*t37.*(t27-t46);
            t192 = L.*g.*m.*t34.*(t26-t45);
            t193 = t63+t68+t92;
            t197 = t71+t89+t96+t99;
            t198 = t56+t70+t79+t88+t98;
            t205 = t50+t67+t86+t90+t97+t105;
            t103 = cos(t102);
            t104 = sin(t102);
            t119 = -t108;
            t121 = q0.*q1.*t10.*t72;
            t122 = q0.*q1.*t10.*t74;
            t123 = t11.*t72.*1.2e+1;
            t125 = t51.*t72;
            t128 = t11.*t74.*1.2e+1;
            t130 = m.*t51.*t74;
            t133 = q1.*t129;
            t136 = dq0.*m.*t109;
            t137 = q1.*t124;
            t138 = m.*t11.*t74.*-2.4e+1;
            t140 = dq0.*m.*t126;
            t141 = q0.*t10.*t31.*t72;
            t142 = m.*t62.*t72;
            t144 = q1.*t10.*t14.*t74;
            t146 = m.*q1.*t10.*t14.*t72;
            t147 = t16.*t124;
            t149 = dq0.*m.*t62.*t74;
            t150 = dq0.*m.*t134;
            t157 = t10.*t14.*t31.*t74;
            t158 = m.*t16.*t129;
            t160 = m.*t10.*t14.*t31.*t72;
            t195 = t62+t66+t91+t95;
            t196 = (dtht1.*m.*t39.*t193)./1.2e+1;
            t200 = t64+t69+t87+t93+t106;
            t201 = (dq1.*m.*t42.*t197)./3.6e+1;
            t202 = (m.*t41.*t198)./2.4e+1;
            t207 = (m.*t38.*t205)./8.0;
            t132 = q0.*t128;
            t135 = -t128;
            t139 = dq0.*m.*t125;
            t151 = q0.*q1.*t10.*t103;
            t152 = q0.*q1.*t10.*t104;
            t153 = t11.*t103.*1.2e+1;
            t154 = t51.*t103;
            t156 = t11.*t104.*1.2e+1;
            t162 = t11.*t103.*-2.4e+1;
            t164 = m.*t51.*t104;
            t168 = t63.*t104;
            t170 = m.*t62.*t103;
            t171 = m.*t63.*t103;
            t173 = q0.*t10.*t31.*t104;
            t174 = q0.*t10.*t16.*t104;
            t175 = q1.*t10.*t14.*t104;
            t176 = q0.*t11.*t104.*-1.2e+1;
            t179 = q0.*q1.*t11.*t103.*-1.2e+1;
            t180 = dq0.*m.*t62.*t104;
            t182 = m.*q1.*t62.*t104;
            t185 = m.*q1.*t10.*t14.*t103;
            t190 = m.*t10.*t14.*t16.*t104;
            t191 = m.*t11.*t16.*t104.*-1.2e+1;
            t199 = (dtht0.*m.*t36.*t195)./6.0;
            t203 = dq1.*t202;
            t204 = dtht1.*t202;
            t206 = (dq0.*m.*t40.*t200)./9.0;
            t208 = dq0.*t207;
            t209 = dtht0.*t207;
            t161 = -t153;
            t163 = q1.*t153;
            t166 = q0.*t156;
            t167 = q1.*t156;
            t169 = dq0.*m.*t154;
            t177 = t16.*t153;
            t178 = m.*q0.*t162;
            t181 = dq0.*m.*t168;
            t186 = dq0.*q1.*t170;
            t187 = -t174;
            t188 = m.*t16.*t156;
            t194 = -t190;
            t172 = q0.*t163;
            t189 = m.*q0.*t177;
            t210 = t111+t122+t129+t135+t145+t156+t163+t173;
            t214 = t119+t123+t124+t132+t133+t157+t161+t167+t175+t176+t179;
            t215 = t120+t125+t126+t134+t141+t151+t155+t162+t168+t177+t187;
            t211 = (m.*t33.*t37.*t210)./1.2e+1;
            t216 = (m.*t34.*t37.*t214)./1.2e+1;
            t220 = (m.*t33.*t39.*t215)./1.2e+1;
            t221 = dq0.*m.*t34.*t37.*t214.*(-1.0./1.2e+1);
            t222 = dtht1.*m.*t34.*t37.*t214.*(-1.0./1.2e+1);
            t225 = t127+t131+t138+t142+t148+t160+t164+t165+t171+t178+t182+t185+t189+t191+t194;
            t212 = dtht0.*t211;
            t213 = dtht1.*t211;
            t217 = dq0.*t216;
            t218 = dtht1.*t216;
            t219 = -t216;
            t223 = dq1.*t220;
            t224 = dtht0.*t220;
            t226 = (t34.*t39.*t225)./1.2e+1;
            M = reshape([(m.*t40.*t200)./9.0,t207,t226,t219,t207,(m.*t36.*t195)./6.0,t220,t211,t226,t220,(m.*t42.*t197)./3.6e+1,t202,t219,t211,t202,(m.*t39.*t193)./1.2e+1],[4,4]);
            if nargout > 1
                et1 = dq1.*t185.*-2.0+dtht1.*t143+dtht1.*t190+dtht1.*t16.*t131+dtht1.*t16.*t164+dtht1.*t16.*t178+dtht1.*m.*q1.*t120+dtht1.*m.*q1.*t162-dq1.*m.*t11.*t60.*7.2e+1+dq1.*m.*t11.*t74.*7.2e+1-dq1.*m.*t11.*t75.*7.2e+1-dq1.*m.*t11.*t104.*7.2e+1+dq1.*m.*t75.*t88+dq1.*m.*t18.*t153+dq1.*m.*t18.*t166+dq1.*m.*t12.*t175+dq1.*m.*t88.*t104+dtht1.*m.*t63.*t72+dtht1.*m.*t18.*t153+dtht1.*m.*t18.*t166-dq1.*m.*q0.*t11.*t72.*7.2e+1+dq1.*m.*q1.*t11.*t73.*7.2e+1+dq1.*m.*q0.*t11.*t103.*7.2e+1-dq1.*m.*q1.*t11.*t103.*7.2e+1+dtht1.*m.*q1.*t62.*t74-dq1.*m.*t11.*t18.*t73.*1.2e+1+dq1.*m.*t14.*t68.*t103-dtht1.*m.*t11.*t18.*t73.*1.2e+1+dtht1.*m.*t14.*t68.*t103+dq1.*m.*t10.*t12.*t14.*t72-dtht1.*m.*t10.*t14.*t16.*t74;
                et2 = dq1.*m.*q0.*q1.*t11.*t104.*-7.2e+1-dtht1.*m.*q0.*q1.*t11.*t104.*2.4e+1-dq1.*m.*q0.*t11.*t16.*t103.*3.6e+1;
                et3 = dtht0.*t142+dtht0.*t178+dtht0.*t182+dtht0.*t189+dq0.*t14.*t130+dq0.*t16.*t170+dtht0.*t14.*t130+dtht0.*m.*q0.*t120+dtht0.*m.*q0.*t155+dq0.*m.*t11.*t60.*4.8e+1-dq0.*m.*t11.*t74.*4.8e+1+dq0.*m.*t11.*t75.*4.8e+1+dq0.*m.*t11.*t104.*4.8e+1+dtht0.*m.*t62.*t73+dq0.*m.*q0.*t11.*t72.*4.8e+1-dq0.*m.*q1.*t11.*t73.*4.8e+1-dq0.*m.*q0.*t11.*t103.*4.8e+1+dq0.*m.*q1.*t11.*t103.*4.8e+1+dq0.*m.*q1.*t66.*t104+dq0.*m.*q1.*t14.*t162+dtht0.*m.*q1.*t62.*t75+dtht0.*m.*q1.*t66.*t104+dtht0.*m.*q1.*t14.*t162-dq0.*m.*t11.*t16.*t75.*2.4e+1-dq0.*m.*t11.*t14.*t104.*2.4e+1-dq0.*m.*t11.*t16.*t104.*2.4e+1+dq0.*m.*t31.*t66.*t74+dq0.*m.*t16.*t66.*t103+dq0.*m.*t16.*t86.*t104-dtht0.*m.*t11.*t14.*t104.*2.4e+1+dtht0.*m.*t31.*t66.*t74+dtht0.*m.*t16.*t66.*t103;
                et4 = dtht0.*m.*t16.*t86.*t104+dq0.*m.*q0.*q1.*t11.*t104.*4.8e+1;
                mt1 = [dq0.*m.*t11.*t34.^3.*(q0.*2.4e+1-t6.*3.0e+1+t24.*6.0-t52.*9.0+t3.*t15.*3.0+t13.*t14).*(-1.0./3.0),m.*t11.*t40.*(dq0.*-4.0+dq0.*t3.*4.0-dq0.*t14.*3.0-dtht0.*t14.*4.0+dtht0.*q0.*t22+dq0.*t3.*t14+dq0.*t6.*t13+dq0.*t6.*t15+dtht0.*t3.*t14+dtht0.*t6.*t15),t36.*t39.*(et3+et4).*(-1.0./1.2e+1)];
                mt2 = [t36.*t37.*(t136+t169+t180+t186+dtht0.*m.*t166+dtht0.*m.*t172+dq0.*q0.*t138+dtht0.*m.*q0.*t111+dtht0.*m.*q0.*t129+dtht0.*m.*q0.*t145-dq0.*m.*t11.*t72.*2.4e+1-dq0.*m.*t11.*t73.*2.4e+1+dq0.*m.*t72.*t86+dtht0.*m.*t72.*t86-dq0.*m.*q1.*t11.*t75.*2.4e+1-dq0.*m.*q1.*t11.*t104.*2.4e+1+dq0.*m.*q1.*t66.*t103+dq0.*m.*q1.*t86.*t104-dtht0.*m.*q0.*t11.*t74.*1.2e+1+dtht0.*m.*q1.*t66.*t103+dtht0.*m.*q1.*t86.*t104-dq0.*m.*t11.*t14.*t103.*1.2e+1+dq0.*m.*t31.*t66.*t72-dtht0.*m.*t11.*t14.*t103.*1.2e+1+dtht0.*m.*t31.*t66.*t72).*(-1.0./1.2e+1),-dtht0.*m.*t11.*t38.*t159,dq0.*m.*t11.*t38.*t159];
                mt3 = [t34.*t39.*(t139+t140+t149+t150+t181+dq0.*m.*t120+dq0.*m.*t155+dq0.*m.*t157+dq0.*m.*t162+dq0.*m.*t175+dq0.*m.*t177+dtht0.*m.*t157+dtht0.*m.*t175+dq0.*q1.*t178+dtht0.*q0.*t158+dtht0.*q1.*t178+dq0.*m.*t16.*t166+dtht0.*m.*t62.*t74+dtht0.*m.*t16.*t166-dq0.*m.*q0.*t11.*t104.*2.4e+1-dtht0.*m.*q0.*t11.*t60.*2.4e+1-dtht0.*m.*q0.*t11.*t75.*2.4e+1-dtht0.*m.*q0.*t11.*t104.*2.4e+1+dtht0.*m.*q1.*t62.*t73+dq0.*m.*t10.*t14.*t16.*t103+dtht0.*m.*t10.*t14.*t16.*t103).*(-1.0./1.2e+1)];
                mt4 = [t34.*t37.*(dq0.*t160+dq0.*t185+dtht0.*t160+dtht0.*t185+dq0.*m.*t111+dq0.*m.*t129+dq0.*m.*t145+dq0.*m.*t156+dq0.*m.*t163+dq0.*m.*q0.*t123+dq0.*m.*q1.*t166+dtht0.*m.*q0.*t123+dtht0.*m.*q0.*t124+dtht0.*m.*q0.*t133+dtht0.*m.*q1.*t166-dq0.*m.*t11.*t74.*1.2e+1-dq0.*m.*q0.*t11.*t103.*1.2e+1-dtht0.*m.*q0.*t11.*t57.*1.2e+1-dtht0.*m.*q0.*t11.*t103.*1.2e+1).*(-1.0./1.2e+1),(t34.*t41.*(et1+et2))./1.2e+1];
                mt5 = [(t33.*t41.*(dq1.*m.*t151.*-2.0+dtht1.*m.*t174+dq1.*m.*t11.*t57.*7.2e+1-dq1.*m.*t11.*t72.*7.2e+1-dq1.*m.*t11.*t73.*7.2e+1+dq1.*m.*t11.*t103.*7.2e+1+dq1.*m.*t18.*t129+dq1.*m.*t73.*t88+dq1.*m.*t12.*t152+dq1.*m.*t18.*t156+dtht1.*m.*t63.*t74+dtht1.*m.*t16.*t126+dtht1.*m.*t18.*t129+dtht1.*m.*t18.*t156+dtht1.*m.*t16.*t162-dq1.*m.*q1.*t11.*t75.*7.2e+1-dq1.*m.*q1.*t11.*t104.*7.2e+1+dq1.*m.*q0.*t68.*t103-dtht1.*m.*q1.*t11.*t60.*2.4e+1-dtht1.*m.*q1.*t11.*t75.*2.4e+1-dtht1.*m.*q1.*t11.*t104.*2.4e+1+dtht1.*m.*q0.*t68.*t103-dq1.*m.*t11.*t16.*t103.*3.6e+1+dq1.*m.*q0.*t10.*t12.*t72-dtht1.*m.*q0.*t10.*t16.*t74))./1.2e+1];
                mt6 = [dq1.*m.*t11.*t37.^3.*(t7.*-1.5e+1+t25.*1.5e+1+t12.*t16+t16.*t23).*(-2.0./3.0),-m.*t11.*t42.*(dq1.*-4.0+dq1.*t4.*4.0+dq1.*t16.*3.0-dq1.*q1.*t7.*2.0-dtht1.*q1.*t7.*3.0+dtht1.*q1.*t12+dq1.*t4.*t16+dtht1.*t4.*t16)];
                mt7 = [(t34.*t39.*(dq1.*m.*t120+dq1.*m.*t125+dq1.*m.*t126+dq1.*m.*t134+dq1.*m.*t155+dq1.*m.*t157+dq1.*m.*t162+dq1.*m.*t168+dq1.*m.*t175+dq1.*m.*t177+dtht1.*m.*t133+dtht1.*m.*t155+dtht1.*m.*t167+dtht1.*m.*t177+dtht1.*m.*t179+dq1.*q1.*t178+dtht1.*m.*q1.*t111+dq1.*m.*t62.*t74+dq1.*m.*t16.*t166+dtht1.*m.*t16.*t166+dtht1.*m.*q0.*q1.*t123-dq1.*m.*q0.*t11.*t104.*2.4e+1-dtht1.*m.*q1.*t11.*t74.*1.2e+1+dq1.*m.*t10.*t14.*t16.*t103-dtht1.*m.*t10.*t14.*t16.*t72+dtht1.*m.*t10.*t14.*t16.*t103))./1.2e+1];
                mt8 = [(t33.*t39.*(dq1.*t130+dq1.*t143+dq1.*t158+dq1.*t188+dtht1.*t158+dtht1.*t188+dq1.*m.*t152+dtht1.*m.*t137+dq1.*m.*q1.*t162+dtht1.*m.*q1.*t123-dq1.*m.*t11.*t60.*2.4e+1-dq1.*m.*t11.*t75.*2.4e+1-dq1.*m.*t11.*t104.*2.4e+1-dtht1.*m.*q1.*t11.*t57.*1.2e+1-dtht1.*m.*q1.*t11.*t103.*1.2e+1+dq1.*m.*q0.*t10.*t31.*t74+dq1.*m.*q0.*t10.*t16.*t103-dtht1.*m.*q0.*t10.*t16.*t72+dtht1.*m.*q0.*t10.*t16.*t103))./1.2e+1,dtht1.*m.*t11.*t41.*t110,-dq1.*m.*t11.*t41.*t110];
                C = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8],4,4);
            end
            if nargout > 2
                G = [t112+k.*q0+L.*g.*m.*t30.*t34.*2.0-L.*g.*m.*t36.*(t26-t45).*2.0;t112-L.*g.*m.*t26.*t33.*2.0-L.*g.*m.*t34.*(t30-t47);t113+k.*q1-L.*g.*m.*t37.*t48.*2.0+L.*g.*m.*t39.*(t27-t46).*2.0;t113-L.*g.*m.*t37.*(t48-sin(t9))];
            end
            if nargout > 3
                t227 = dq0.*t226;
                t228 = dq1.*t226;
                t229 = t196+t203+t212+t221;
                t231 = t199+t208+t213+t223;
                t230 = (dtht1.*t229)./2.0;
                t232 = (dtht0.*t231)./2.0;
                t233 = t201+t204+t224+t227;
                t235 = t206+t209+t222+t228;
                t234 = (dq1.*t233)./2.0;
                t236 = (dq0.*t235)./2.0;
                E = t61+t65+t76+t77+t114+t115+t118+t184+t192+t230+t232+t234+t236;
            end
            if nargout > 4
                PE = t61+t65+t76+t77+t114+t115+t118+t184+t192;
            end
            if nargout > 5
                KE = t230+t232+t234+t236;
            end
        end

        function [M,C,G,E,PE,KE] = getDynsSpring(Robot,q,dq,p)
            %helper function to numerically calculate dynamics matrices M (inertia), C
            %(coriolis) and G (gravity) given state vector q and its derivatives along
            %with hip pos'n vector p
            %this version includes spring potential energy in the PE
            %calculation (PE = gravity + spring) which may cause issues
            %witht eh group action trying to project a spring stiffness
            %that doesn't exist (make sure the group action offset is only
            %applied to the base angles, not curavtures)

            q0 = q(1); tht0 = q(2); q1 = q(3); tht1 = q(4);
            dq0 = dq(1); dtht0 = dq(2); dq1 = dq(3); dtht1 = dq(4);
            footX = p(1); footY = p(2);
            [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = Robot.getVars();

            t2 = cos(phi);
            t3 = cos(q0);
            t4 = cos(q1);
            t5 = sin(phi);
            t6 = sin(q0);
            t7 = sin(q1);
            t8 = phi+tht0;
            t9 = phi+tht1;
            t10 = D.^2;
            t11 = L.^2;
            t12 = q1.*2.0;
            t13 = q0.*4.0;
            t14 = q0.^2;
            t15 = q0.^3;
            t16 = q1.^2;
            t18 = q1.^3;
            t19 = q0.^5;
            t21 = q1.^5;
            t32 = -q1;
            t34 = 1.0./q0;
            t36 = 1.0./q1;
            t44 = -tht0;
            t45 = -tht1;
            t17 = t14.^2;
            t20 = t16.^2;
            t22 = t6.*3.0;
            t23 = t7.*3.0;
            t24 = q0.*t3;
            t25 = q1.*t4;
            t26 = cos(t8);
            t27 = cos(t9);
            t28 = q0+t8;
            t29 = q1+t9;
            t30 = sin(t8);
            t31 = sin(t9);
            t33 = -t13;
            t35 = 1.0./t14;
            t37 = 1.0./t15;
            t38 = 1.0./t16;
            t40 = 1.0./t18;
            t41 = 1.0./t19;
            t43 = 1.0./t21;
            t51 = t11.*8.0;
            t52 = t11.*2.4e+1;
            t53 = t6.*t14;
            t54 = t45+tht0;
            t59 = q1+t44+tht1;
            t61 = footY.*g.*m.*t2.*2.0;
            t64 = q0.*t11.*3.6e+1;
            t65 = footX.*g.*m.*t5.*2.0;
            t66 = t10.*t15;
            t68 = t10.*t18;
            t69 = t10.*t19;
            t71 = t10.*t21;
            t72 = q1.*t11.*7.2e+1;
            t77 = (k.*t14)./2.0;
            t78 = (k.*t16)./2.0;
            t81 = t6.*t11.*1.2e+1;
            t83 = t6.*t11.*3.6e+1;
            t87 = t11.*t14.*1.2e+1;
            t88 = t11.*t15.*1.2e+1;
            t89 = t11.*t16.*1.2e+1;
            t90 = t11.*t18.*1.2e+1;
            t91 = t3.*t11.*-8.0;
            t92 = t4.*t11.*-2.4e+1;
            t94 = t7.*t11.*-2.4e+1;
            t96 = t7.*t11.*1.44e+2;
            t99 = q0.*t6.*t11.*-8.0;
            t39 = 1.0./t17;
            t42 = 1.0./t20;
            t46 = cos(t28);
            t47 = cos(t29);
            t48 = sin(t28);
            t49 = sin(t29);
            t50 = -t23;
            t55 = -t26;
            t56 = -t27;
            t57 = cos(t54);
            t58 = q0+t54;
            t60 = sin(t54);
            t62 = q0.*t52;
            t63 = q1.*t52;
            t67 = t10.*t17;
            t70 = t10.*t20;
            t74 = cos(t59);
            t76 = sin(t59);
            t79 = t3.*t51;
            t80 = t4.*t52;
            t82 = t7.*t52;
            t84 = t11.*t24.*1.2e+1;
            t85 = q0.*t6.*t51;
            t93 = -t81;
            t95 = -t83;
            t98 = t11.*t25.*7.2e+1;
            t100 = q1.*t94;
            t101 = -t96;
            t103 = t11.*t53.*1.8e+1;
            t105 = L.*g.*m.*t31.*t36;
            t108 = t14.*t91;
            t110 = L.*g.*m.*t30.*t34.*2.0;
            t178 = t22+t24+t33+t53;
            t73 = cos(t58);
            t75 = sin(t58);
            t86 = t7.*t63;
            t97 = -t84;
            t102 = t14.*t79;
            t104 = t32+t58;
            t109 = -t103;
            t111 = t11.*t57.*1.2e+1;
            t112 = t52.*t57;
            t113 = t12+t25+t50;
            t114 = t11.*t60.*1.2e+1;
            t115 = t52.*t60;
            t116 = L.*g.*m.*t34.*t46;
            t117 = L.*g.*m.*t34.*t48;
            t118 = t46+t55;
            t119 = t47+t56;
            t120 = -t110;
            t121 = -t105;
            t122 = t11.*t57.*-2.4e+1;
            t126 = t11.*t74.*1.2e+1;
            t129 = t52.*t74;
            t131 = t11.*t76.*1.2e+1;
            t134 = t52.*t76;
            t140 = q1.*t11.*t57.*-1.2e+1;
            t142 = q1.*t11.*t60.*-1.2e+1;
            t157 = q1.*t11.*t74.*-1.2e+1;
            t159 = dq0.*m.*t11.*t74.*-2.4e+1;
            t171 = dtht1.*m.*t57.*t89;
            t172 = dtht1.*m.*t60.*t89;
            t181 = dtht1.*m.*t11.*t16.*t57.*-1.2e+1;
            t195 = -L.*g.*m.*t35.*(t26-t46);
            t196 = -L.*g.*m.*t38.*(t27-t47);
            t199 = L.*g.*m.*t35.*(t26-t46);
            t200 = L.*g.*m.*t38.*(t27-t47);
            t201 = t63+t68+t94;
            t204 = t52+t70+t89+t92+t100;
            t207 = t71+t72+t90+t98+t101;
            t211 = t51+t67+t87+t91+t99+t108;
            t106 = cos(t104);
            t107 = sin(t104);
            t123 = q0.*q1.*t10.*t73;
            t124 = q0.*q1.*t10.*t75;
            t125 = t11.*t73.*1.2e+1;
            t127 = q1.*t111;
            t128 = t52.*t73;
            t130 = t11.*t75.*1.2e+1;
            t132 = q1.*t114;
            t133 = t52.*t75;
            t137 = q1.*t131;
            t139 = -t126;
            t143 = t11.*t75.*-2.4e+1;
            t144 = dq0.*m.*t112;
            t146 = q1.*t126;
            t147 = t62.*t73;
            t148 = q0.*t11.*t75.*-1.2e+1;
            t149 = q1.*t11.*t75.*-1.2e+1;
            t151 = dq0.*m.*t129;
            t153 = q0.*t10.*t32.*t75;
            t155 = q1.*t10.*t14.*t73;
            t156 = q1.*t10.*t14.*t75;
            t158 = dq0.*m.*t11.*t73.*-2.4e+1;
            t161 = dq0.*m.*t62.*t75;
            t173 = dtht0.*m.*q0.*t140;
            t174 = t10.*t14.*t32.*t75;
            t176 = dtht0.*m.*q0.*t142;
            t177 = dtht1.*m.*t73.*t89;
            t179 = dtht1.*m.*t75.*t89;
            t188 = dtht1.*m.*t11.*t16.*t75.*-1.2e+1;
            t202 = t62+t66+t93+t97;
            t203 = (dtht1.*m.*t40.*t201)./1.2e+1;
            t206 = t64+t69+t88+t95+t109;
            t208 = (m.*t42.*t204)./2.4e+1;
            t213 = (dq1.*m.*t43.*t207)./3.6e+1;
            t214 = (m.*t39.*t211)./8.0;
            t135 = q0.*t130;
            t136 = q1.*t130;
            t138 = -t125;
            t141 = -t130;
            t145 = q1.*t125;
            t150 = dq0.*m.*t128;
            t160 = dtht0.*m.*q0.*t127;
            t162 = dtht0.*m.*q0.*t132;
            t163 = q0.*q1.*t10.*t106;
            t165 = q0.*q1.*t10.*t107;
            t167 = t11.*t106.*1.2e+1;
            t168 = t52.*t106;
            t169 = t11.*t107.*1.2e+1;
            t170 = t52.*t107;
            t175 = dq0.*m.*q0.*t143;
            t180 = t11.*t106.*-2.4e+1;
            t183 = t62.*t106;
            t187 = q0.*t10.*t32.*t106;
            t190 = q1.*t10.*t14.*t106;
            t191 = q0.*t10.*t16.*t107;
            t192 = q1.*t10.*t14.*t107;
            t194 = dq0.*m.*t62.*t107;
            t197 = t10.*t14.*t16.*t107;
            t198 = t10.*t14.*t32.*t106;
            t205 = (dtht0.*m.*t37.*t202)./6.0;
            t209 = dq1.*t208;
            t210 = dtht1.*t208;
            t212 = (dq0.*m.*t41.*t206)./9.0;
            t215 = dq0.*t214;
            t216 = dtht0.*t214;
            t152 = q0.*t145;
            t154 = q1.*t135;
            t182 = q1.*t167;
            t184 = q0.*t169;
            t185 = q1.*t169;
            t186 = dq0.*m.*t168;
            t193 = q0.*t180;
            t217 = t114+t131+t140+t141+t145+t153+t165+t169;
            t164 = dtht0.*m.*t152;
            t166 = dtht0.*m.*t154;
            t189 = q1.*t184;
            t218 = (m.*t34.*t38.*t217)./1.2e+1;
            t222 = dtht0.*m.*t34.*t38.*t217.*(-1.0./1.2e+1);
            t223 = dtht1.*m.*t34.*t38.*t217.*(-1.0./1.2e+1);
            t224 = t122+t123+t128+t129+t136+t137+t142+t180+t185+t187+t191;
            t225 = t111+t132+t138+t139+t148+t149+t152+t167+t174+t184+t192;
            t219 = dtht0.*t218;
            t220 = dtht1.*t218;
            t221 = -t218;
            t226 = (m.*t34.*t40.*t224)./1.2e+1;
            t229 = (m.*t35.*t38.*t225)./1.2e+1;
            t233 = dq1.*m.*t34.*t40.*t224.*(-1.0./1.2e+1);
            t234 = dtht0.*m.*t34.*t40.*t224.*(-1.0./1.2e+1);
            t236 = dq0.*m.*t35.*t38.*t225.*(-1.0./1.2e+1);
            t237 = dtht1.*m.*t35.*t38.*t225.*(-1.0./1.2e+1);
            t238 = t115+t134+t140+t143+t145+t147+t154+t155+t157+t170+t182+t189+t193+t197+t198;
            t227 = dq1.*t226;
            t228 = dtht0.*t226;
            t230 = -t226;
            t231 = dq0.*t229;
            t232 = dtht1.*t229;
            t235 = -t229;
            t239 = (m.*t35.*t40.*t238)./1.2e+1;
            t243 = dq0.*m.*t35.*t40.*t238.*(-1.0./1.2e+1);
            t244 = dq1.*m.*t35.*t40.*t238.*(-1.0./1.2e+1);
            t245 = t203+t209+t222+t236;
            t247 = t205+t215+t223+t233;
            t240 = dq0.*t239;
            t241 = dq1.*t239;
            t242 = -t239;
            M = reshape([(m.*t41.*t206)./9.0,t214,t242,t235,t214,(m.*t37.*t202)./6.0,t230,t221,t242,t230,(m.*t43.*t207)./3.6e+1,t208,t235,t221,t208,(m.*t40.*t201)./1.2e+1],[4,4]);
            if nargout > 1
                et1 = t172+t188+q0.*t177-dq1.*m.*t190.*2.0+dtht1.*m.*t197+dq1.*m.*q1.*t122+dq1.*m.*t11.*t60.*7.2e+1-dq1.*m.*t11.*t75.*7.2e+1+dq1.*m.*t11.*t76.*7.2e+1+dq1.*m.*t11.*t107.*7.2e+1+dq1.*m.*t63.*t73+dq1.*m.*t12.*t192+dtht1.*m.*t57.*t63+dtht1.*m.*t63.*t106+dtht1.*m.*q0.*q1.*t143+dq1.*m.*q0.*t11.*t73.*7.2e+1-dq1.*m.*q1.*t11.*t74.*4.8e+1-dq1.*m.*q0.*t11.*t106.*7.2e+1+dq1.*m.*q1.*t11.*t106.*4.8e+1+dq1.*m.*q1.*t62.*t75+dq1.*m.*q0.*t89.*t106-dtht1.*m.*q1.*t11.*t73.*2.4e+1-dtht1.*m.*q1.*t11.*t74.*2.4e+1+dtht1.*m.*q1.*t62.*t107+dtht1.*m.*q0.*t89.*t106-dq1.*m.*t11.*t16.*t76.*1.2e+1-dq1.*m.*t11.*t16.*t107.*1.2e+1+dq1.*m.*t14.*t68.*t106-dtht1.*m.*t11.*t16.*t76.*1.2e+1-dtht1.*m.*t11.*t16.*t107.*1.2e+1;
                et2 = dtht1.*m.*t14.*t68.*t106+dq1.*m.*t10.*t12.*t14.*t73-dtht1.*m.*t10.*t14.*t16.*t75+dq1.*m.*q0.*q1.*t11.*t107.*4.8e+1;
                et3 = t166+t176+q1.*t159+q1.*t161+q1.*t194+dtht0.*m.*t147+dtht0.*m.*t189+dtht0.*m.*t193+dq0.*m.*q1.*t122+dtht0.*m.*q0.*t122+dtht0.*m.*q0.*t137+dq0.*m.*t11.*t60.*4.8e+1-dq0.*m.*t11.*t75.*4.8e+1+dq0.*m.*t11.*t76.*4.8e+1+dq0.*m.*t11.*t107.*4.8e+1+dq0.*m.*t63.*t73+dq0.*m.*t14.*t133+dq0.*m.*t63.*t106+dtht0.*m.*t62.*t74+dtht0.*m.*t14.*t133+dq0.*m.*q0.*t11.*t73.*4.8e+1-dq0.*m.*q0.*t11.*t106.*4.8e+1+dq0.*m.*q1.*t66.*t75+dtht0.*m.*q1.*t66.*t75-dq0.*m.*t11.*t14.*t107.*2.4e+1-dq0.*m.*t16.*t66.*t106+dq0.*m.*t32.*t66.*t107-dtht0.*m.*t11.*t14.*t107.*2.4e+1-dtht0.*m.*t16.*t66.*t106+dtht0.*m.*t32.*t66.*t107-dq0.*m.*q1.*t11.*t14.*t73.*1.2e+1-dq0.*m.*q1.*t11.*t14.*t106.*1.2e+1-dtht0.*m.*q1.*t11.*t14.*t73.*1.2e+1;
                et4 = dtht0.*m.*q1.*t11.*t14.*t106.*-1.2e+1;
                mt1 = [dq0.*m.*t11.*t35.^3.*(q0.*2.4e+1-t6.*3.0e+1+t24.*6.0-t53.*9.0+t3.*t15.*3.0+t13.*t14).*(-1.0./3.0),m.*t11.*t41.*(dq0.*-4.0+dq0.*t3.*4.0-dq0.*t14.*3.0-dtht0.*t14.*4.0+dtht0.*q0.*t22+dq0.*t3.*t14+dq0.*t6.*t13+dq0.*t6.*t15+dtht0.*t3.*t14+dtht0.*t6.*t15),(t37.*t40.*(et3+et4))./1.2e+1];
                mt2 = [(t37.*t38.*(t144+t158+t159+t164+t173+t175+t186+t194+dtht0.*m.*t148+dtht0.*m.*t184+dq0.*m.*q1.*t143+dq0.*m.*q1.*t147+dtht0.*m.*q0.*t114+dtht0.*m.*q0.*t131+dq0.*m.*t60.*t63+dq0.*m.*t73.*t87+dtht0.*m.*t73.*t87+dq0.*m.*q1.*t66.*t73+dq0.*m.*q1.*t75.*t87+dtht0.*m.*q1.*t66.*t73+dtht0.*m.*q1.*t75.*t87-dq0.*m.*t11.*t14.*t106.*1.2e+1+dq0.*m.*t32.*t66.*t106-dtht0.*m.*t11.*t14.*t106.*1.2e+1+dtht0.*m.*t32.*t66.*t106))./1.2e+1,-dtht0.*m.*t11.*t39.*t178,dq0.*m.*t11.*t39.*t178];
                mt3 = [t35.*t40.*(t144+t158+t159+t164+t173+t175+t186+t194+dq0.*m.*t132+dq0.*m.*t149+dq0.*m.*t152+dq0.*m.*t174+dq0.*m.*t192+dtht0.*m.*t174+dtht0.*m.*t192+dq0.*m.*q0.*t182+dtht0.*m.*q0.*t143+dtht0.*m.*q0.*t157+dtht0.*m.*q0.*t182+dtht0.*m.*t60.*t62+dtht0.*m.*t62.*t76+dtht0.*m.*t62.*t107-dq0.*m.*q1.*t11.*t76.*1.2e+1-dq0.*m.*q1.*t11.*t107.*1.2e+1+dq0.*m.*t10.*t14.*t16.*t106+dtht0.*m.*t10.*t14.*t16.*t106).*(-1.0./1.2e+1)];
                mt4 = [(t35.*t38.*(t166+t176+dq0.*m.*t114+dq0.*m.*t131+dq0.*m.*t140+dq0.*m.*t145+dq0.*m.*t154+dq0.*m.*t155+dq0.*m.*t169+dq0.*m.*t198+dtht0.*m.*t155+dtht0.*m.*t198+dq0.*m.*q0.*t125+dtht0.*m.*q0.*t125+dtht0.*m.*q0.*t126-dq0.*m.*t11.*t75.*1.2e+1-dq0.*m.*q0.*t11.*t106.*1.2e+1-dtht0.*m.*q0.*t11.*t57.*1.2e+1-dtht0.*m.*q0.*t11.*t106.*1.2e+1))./1.2e+1,(t35.*t42.*(et1+et2))./1.2e+1];
                mt5 = [(t34.*t42.*(t177+t181-dq1.*m.*t163.*2.0+dtht1.*m.*t191+dtht1.*m.*q1.*t143-dq1.*m.*t11.*t57.*7.2e+1+dq1.*m.*t11.*t73.*7.2e+1+dq1.*m.*t11.*t74.*7.2e+1-dq1.*m.*t11.*t106.*7.2e+1+dq1.*m.*t63.*t75+dq1.*m.*t12.*t165+dq1.*m.*t89.*t106+dtht1.*m.*t60.*t63+dtht1.*m.*t63.*t76+dtht1.*m.*t63.*t107+dtht1.*m.*t89.*t106-dq1.*m.*q1.*t11.*t60.*2.4e+1+dq1.*m.*q1.*t11.*t76.*4.8e+1+dq1.*m.*q1.*t11.*t107.*4.8e+1+dq1.*m.*q0.*t68.*t106+dtht1.*m.*q0.*t68.*t106-dq1.*m.*t11.*t16.*t74.*1.2e+1-dtht1.*m.*t11.*t16.*t74.*1.2e+1+dq1.*m.*q0.*t10.*t12.*t73-dtht1.*m.*q0.*t10.*t16.*t75))./1.2e+1,dq1.*m.*t11.*t38.^3.*(q1.*1.2e+1-t7.*3.0e+1+t18+t25.*1.8e+1+t16.*t23).*(-1.0./3.0)];
                mt6 = [-m.*t11.*t43.*(dq1.*4.0-dq1.*t4.*4.0+dq1.*t16-dq1.*q1.*t7.*4.0-dtht1.*q1.*t7.*3.0+dtht1.*q1.*t12+dq1.*t4.*t16+dtht1.*t4.*t16),t35.*t40.*(t177+t181+q0.*t179+dq1.*m.*t122+dq1.*m.*t128+dq1.*m.*t129+dq1.*m.*t136+dq1.*m.*t137+dq1.*m.*t142+dq1.*m.*t156+dq1.*m.*t180+dq1.*m.*t185+dtht1.*m.*t132+dtht1.*m.*t137+dtht1.*m.*t149+dtht1.*m.*t152+dtht1.*m.*t185+dq1.*m.*t62.*t75-dq1.*m.*q0.*t11.*t107.*2.4e+1-dq1.*m.*t10.*t14.*t16.*t106+dq1.*m.*t10.*t14.*t32.*t107+dtht1.*m.*t10.*t14.*t16.*t73-dtht1.*m.*t10.*t14.*t16.*t106-dq1.*m.*q0.*q1.*t11.*t73.*1.2e+1-dq1.*m.*q0.*q1.*t11.*t106.*1.2e+1-dtht1.*m.*q0.*q1.*t11.*t106.*1.2e+1).*(-1.0./1.2e+1)];
                mt7 = [(t34.*t40.*(t172+t188+dq1.*m.*t115+dq1.*m.*t134+dq1.*m.*t140+dq1.*m.*t143+dq1.*m.*t145+dq1.*m.*t153+dq1.*m.*t157+dq1.*m.*t165+dq1.*m.*t170+dq1.*m.*t182+dtht1.*m.*t127+dtht1.*m.*t157+dtht1.*m.*t182-dtht1.*m.*q1.*t11.*t73.*1.2e+1+dq1.*m.*q0.*t10.*t16.*t106-dtht1.*m.*q0.*t10.*t16.*t73+dtht1.*m.*q0.*t10.*t16.*t106))./1.2e+1,dtht1.*m.*t11.*t42.*t113,-dq1.*m.*t11.*t42.*t113];
                C = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],4,4);
            end
            if nargout > 2
                G = [t116+k.*q0+L.*g.*m.*t30.*t35.*2.0-L.*g.*m.*t37.*(t26-t46).*2.0;t116-L.*g.*m.*t26.*t34.*2.0-L.*g.*m.*t35.*(t30-t48);k.*q1+L.*g.*m.*t31.*t38+L.*g.*m.*t38.*t49-L.*g.*m.*t40.*(t27-t47).*2.0;L.*g.*m.*t36.*t56-L.*g.*m.*t38.*(t31-t49)];
            end
            if nargout > 3
                t246 = (dtht1.*t245)./2.0;
                t248 = (dtht0.*t247)./2.0;
                t249 = t210+t213+t234+t243;
                t251 = t212+t216+t237+t244;
                t250 = (dq1.*t249)./2.0;
                t252 = (dq0.*t251)./2.0;
                E = t61+t65+t77+t78+t117+t120+t121+t199+t200+t246+t248+t250+t252;
            end
            if nargout > 4
                PE = t61+t65+t77+t78+t117+t120+t121+t199+t200;
            end
            if nargout > 5
                KE = t246+t248+t250+t252;
            end
        end

        function [D,L,g,m,phi,dt,beta,ke,kd,kp,k,nPoints,alpha,epsilon,tht0des] = getVars(Robot)
            %helper function; returns state vars for computations elsewhere
            %D - arm width
            %L - arm length
            %g - |gravity|
            %m - linkage mass (by section?)
            %phi - g field angle from negative vertical, ccw

            D = .1;
            L = 1;
            g = 9.81;
            m = 2;
            phi = 0;
            dt = .01;
            beta = 0; %.1 default
            ke = 3;
            kd = 1.25;
            kp = 1.1;
            k = 1.72; %see sheet of spring constants
            nPoints = 15;
            alpha = .9;
            epsilon = .1;
            tht0des = pi/4;
        end

        function drawStates(Robot,y,t,ye,te)
            %should find a better way of sanitizing this
            q = y(:,1:4); dq = y(:,5:8);
            figure("Name","State Space Curves","Position",[0 250 500 500]);
            subplot(2,2,1); hold on; title("Member Curvatures");
            plot(q(:,1),q(:,3));
            xlabel("q_{st}"); ylabel("q_{sw}");
            subplot(2,2,2); hold on; title("Member Curvature Derivatives");
            plot(dq(:,1),dq(:,3));
            xlabel("dq_{st}"); ylabel("dq_{sw}");
            subplot(2,2,3); hold on; title("Base Angles");
            plot(q(:,2),q(:,4));
            xlabel("\theta_{st}"); ylabel("\theta_{sw}");
            subplot(2,2,4); hold on; title("Base Angle Derivatives");
            plot(dq(:,2),dq(:,4));
            xlabel("d\theta_{st}"); ylabel("d\theta_{sw}");

            figure("Name", "Time Evolution of System Parameters","Position",[0 0 1000 500]);
            subplot(2,4,1); hold on; title("Stance Leg Curvature");
            plot(t,q(:,1));
            xlabel("t"); ylabel("q_{st}");
            subplot(2,4,5); hold on; title("Stance Leg Curvature Derivative");
            plot(t,dq(:,1));
            xlabel("t"); ylabel("dq_{st}");
            subplot(2,4,2); hold on; title("Stance Leg Base Angle");
            plot(t,q(:,2));
            xlabel("t"); ylabel("\theta_{st}");
            subplot(2,4,6); hold on; title("Stance Leg Base Angle Derivative");
            plot(t,dq(:,2));
            xlabel("t"); ylabel("d\theta_{st}");
            subplot(2,4,3); hold on; title("Swing Leg Curvature");
            plot(t,q(:,3));
            xlabel("t"); ylabel("q_{sw}");
            subplot(2,4,7); hold on; title("Swing Leg Curvature Derivative");
            plot(t,dq(:,3));
            xlabel("t"); ylabel("dq_{sw}");
            subplot(2,4,4); hold on; title("Swing Leg Base Angle");
            plot(t,q(:,4));
            xlabel("t"); ylabel("\theta_{st}");
            subplot(2,4,8); hold on; title("Swing Leg Base Angle Derivative");
            plot(t,dq(:,4));
            xlabel("t"); ylabel("d\theta_{st}");
        end

        function drawHeight(Robot,y,t)
            figure("Name","Swing Foot Height")
            h = [];
            for i = 1:size(y,1)
                h = [h Robot.getHeight(y(i,1:4))];
            end
            plot(t,h); hold on;
            xlabel("Time"); ylabel("Swing Foot Height");
            title("Swing Foot Height")
            plot([0 t(end)],[0 0])
        end

    end
end

