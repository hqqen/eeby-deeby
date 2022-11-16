classdef World
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        groundPointsX
        groundPointsY
        ground
    end
    
    methods

        function h = getHeight(World,x,y)
            %unambiguous height-above-ground function for swing foot
            h = y - interpGrnd(World,x);
        end
        
%         function isImpact = findImpact(Robot)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             is = obj.Property1 + inputArg;
%         end

        function phi = getSlope(World,x)
            %get the (piecewise linear) slope at a point; found by taking
            %the arctan of the slope of the pwise linear ground surface at
            %the given point x
            leftInd = find(World.groundPointsX < x,1);
            rightInd = find(World.groundPointsX > x,1);
            phi = atan2(World.groundPointsY(rightInd) - World.groundPointsY(leftInd),...
                        World.groundPointsX(rightInd) - World.groundPointsX(leftInd));

        end

        function h = interpGrnd(World,x)
            %given two points along the ground surface check if they
            %correspond to a point above a ground vertex, then interpoalte
            %if they don't
            %returns gound height at point x
            points = find(World.groundPointsX == x);
            if points ~= []
                %our given point is on a given world point,return that
                %hiehgt difference w no interpolation
                h = y - World.groundPointsY(points);
            else
                leftInd = find(World.groundPointsX < x,1);
                rightInd = find(World.groundPointsX > x,1);
                h = linterp([World.groundPointsX(leftInd), World.groundPointsX(rightInd)],...
                            [World.groundPointsY(leftInd), World.groundPointsY(rightInd)],...
                            x);
            end
        end

        function makeWorld(World)
        %verify the coordinate vectors uesd to build the world are properly
        %structured
            if numel(World.groundPointsX) ~= numel(World.groundPointsY)
                fprintf("Invalid World Coordinate Vectors - Size Mismatch, Probably Missing Entries\n")
                return
            elseif (size(World.groundPointsX,1) ~= 1 && size(World.groundPointsX,2) ~= 1) || (size(World.groundPointsY,1) ~= 1 && size(World.groundPointsY,2) ~= 1) 
                fprintf("Invalid World Coordinate Vectors - Improper Shape, Need 1xN or Nx1 Vector\n")
                return
            end
            World.groundPointsX = reshape(World.groundPointsX,1,numel(World.groundPointsX));
            World.groundPointsY = reshape(World.groundPointsY,1,numel(World.groundPointsY));
            World.ground = sort([World.groundPointsX;World.groundPointsY],2);
            fprintf("World Succesfully Built!\n");
        end
    end
end

