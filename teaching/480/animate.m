function animate(x, a1, a2)
    %given an inegrated state vector x and link lengths a1, a2
    %plot an animted version of the double pendulum

    tht1 = x(:,1); tht2 = x(:,2);
    
    
    figure(); 

    for i = 1:max(size(x))

    p1 = [...
        a1*sin(tht1(i));
        -a1*cos(tht1(i))];

    p2 = p1 + [...
        a2*sin(tht2(i));
        -a2*cos(tht2(i))];

        plot([0, p1(1)],[0, p1(2)], 'b-'); hold on;
        plot([p1(1), p2(1)], [p1(2), p2(2)], 'r-');
        xlim([-a1 - a2, a1 + a2])
        ylim([-a1 - a2, a1 + a2])
        title("Double Pendulum Animation"); hold off; drawnow
    end

end