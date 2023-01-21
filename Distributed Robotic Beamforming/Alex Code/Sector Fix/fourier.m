function [s] = squareFourier(n,P,D,A,x)
    % generate fourier series of square wave up to the nth harmonic
    % INPUTS:
    % n - desired harmonic
    % P - distance ocovered in one period
    % D - wave duty cycle (0<D<1)
    % A - square wave amplitude
    % x - desired points to sample series at
    % OUTPUTS:
    % s - alue of fourier series at sampled points
    % plot of fourier series over one period

    s = A*D; % a0 - this is what will be returned
    xPlot = 0:.01:P; % these two are for plotting and will not be returned
    sPlot = s;
    
    
    for N = 1:n % sum across harmonics
        an = (A/(N*pi))*sin(2*pi*N*D);
        bn = (2*A/(N*pi))*sin(pi*N*D)^2;
        s = s + an.*cos(2.*pi.*N.*x./P) + bn.*sin(2.*pi.*N.*x./P);
        sPlot = sPlot + an.*cos(2.*pi.*N.*xPlot./P) + bn.*sin(2.*pi.*N.*xPlot./P);
    end
    
    figure(90210)
    plot(xPlot,sPlot)
    drawnow

end