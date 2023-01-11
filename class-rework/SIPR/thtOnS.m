function [isOutOfRange,isTerminal,direction] = thtOnS(t,y)
%helper function to stop integration and keep tht on S1 (0,2pi)
isTerminal = 1; direction = 0; isOutOfRange = 0;
if y(2) > 2*pi | y(2) < 0
    isOutOfRange = 1;
end


end