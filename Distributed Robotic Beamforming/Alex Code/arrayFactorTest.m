x = [0, 0];
y = [0, 2*pi];
gain = [1, 1];
a = [1, 1];
alpha = [0, 0];
tht = 0:.01:2*pi;
AF = [tht; zeros(1, size(tht,2))];
for d = 0:pi/100:8*pi
    gain(2) = d;
    for j = 1:max(size(tht))
        for k = 1:max(size(x))

            AF(2,j) = AF(2,j) + (a(k) * exp(i * (alpha(k)  + gain(k)*x(k)*cos(tht(j)) + gain(k)*y(k)*sin(tht(j)))));

        end
    end
    polarplot(AF(1,:), abs(AF(2,:)))
    drawnow
end