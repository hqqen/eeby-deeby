%Fikadu Dagefu, Updated Feb 10 2020   
%Phase front variation simulations
%The scene considered here is a section of Los angeles 
clc; clear;

%% Constants 
lambda = 7.5; 
N = 1;                                    % Number of elements
d = lambda/2;                             % antenna element spacing 
nn = 1;
k = (2.*pi)./lambda; 


%% load simulated FDTD data
grid = load('LA_FDTDin.txt'); 
nx = 1136; %493; %722;                                     %number of grid points in x
ny = 1056; %709; %1019;                                    %number of grid points in y
xgrid = grid(1:nx); ygrid = grid(nx+2:nx+ny+1);            %discretization steps (non-uniform) 
xv = cumsum(xgrid(1:end));  yv = cumsum(ygrid);            %vectors in x and y
xv = xv - round(max(xv)./2) + 30;  %9                         %Adjust the grids to put the source at the correct location (x=??, y=??)
yv = yv - 220; %143                                            %Adjust the grids to put the source at the correct location (x=??, y=??)
[xmat,ymat] = meshgrid(yv,xv);                             %x,y coordinates of receiver locations (signal samples) 
rmat = (xmat.^2 + ymat.^2).^0.5;                           %Range matrix

% save('x_sensor.mat','xv')
% save('y_sensor.mat','yv')

pmat_fs = angle(exp(-1i.*k.*rmat));                         %free space phase map (using the same sampling points as FDTD) in radians

S4g = load('Sensor_3.SEN');  
E4g = S4g(:,8).*exp(1i*S4g(:,9)); 
Emat4g = reshape(E4g,ny,nx).'; Emat4g = Emat4g./max(max((Emat4g))); 
pmat_fdtd = (angle(Emat4g));                                 %phase map from FDTD simulations (in radians)
amat_fdtd = (abs(Emat4g));                                   %amplitude map from FDTD simulations
a = reshape(amat_fdtd,1,1136*1056);
% plot(a)
% save('amat_sensor_7.mat','amat_fdtd')

%% In order to reduce the data size, you can save the pmat_fdtd data 
%save pmat_fdtd pmat_fdtd
%load pmat_fdtd 

% % % %Sanity check, plot the phase front 
figure(1)
pcolor(xv,yv,pmat_fdtd')
shading flat
% 
%% For a given range, compute the phase distribution
%rcr = 100;                                             %Range of interest 
%pe_spatial_deg = 2;%0.3;                                %inherent phase error we can tolerate (due to the fact that we are not sampling on a circle for each range) 
delr = lambda/5; %pe_spatial_deg*(pi/180)/k;  

pmat_fdtd_current = zeros(nx,ny); 
rcrv = 1:delr:max(max(rmat));                               %Range vector    
nf = length(rcrv);                                       %number of frames 
Fw_lw(nf) = struct('cdata',[],'colormap',[]);
W2 = VideoWriter('Phase_front.avi'); W2.FrameRate = 3; open(W2);

for ii = 1:length(rcrv)

    rcr = rcrv(ii);                                      %current range 
    rcs = rcr - delr;                                    %choose the points that are the "same" range (or are on a circle)  
    rcl = rcr + delr;

    idx = find(rmat>=rcs & rmat<=rcl);                   %corresponding indices
    [idxc,idyc] = find(rmat>=rcs & rmat<=rcl);
    
    for jj = 1:length(idxc)
    
        pmat_fdtd_current(idxc(jj),idyc(jj)) = pmat_fdtd(idxc(jj),idyc(jj));
        
    end
    
    figure(2)
    pcolor(yv,xv,pmat_fdtd_current);  shading flat; colorbar;     
    title('Example')
    Fw_lw(ii) = getframe;   
    writeVideo(W2,getframe);
    
end

close(W2);

fps = 10; n = 10; 
movie(Fw_lw,n,fps)

