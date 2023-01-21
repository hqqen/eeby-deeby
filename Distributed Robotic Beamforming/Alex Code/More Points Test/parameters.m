%========================================================================= 
clc
clear all
close all
rng default

%==========================================================================
% Parameters of desired beam
f = 40e6;                                     % Frequency in Hz                                                                        
N_d = 5;                                       % Number of elements of array  
lambda = 3e8/f;                               % wave length
d = lambda/4;                                 % Inter-element spacing
k = (2*pi)./lambda;                           % wave number
% channel parameters
gamma = 1;                                   
mu = 2;
nu = 1;  
rho = 1.5*d; 

                                      
% desired beam (array factor) parameters
%==========================================================================
Am_d = zeros(1,N_d);                                    % element amplitude
rm_d = zeros(2,N_d);                                    % element locations
alpha_d = 45*ones(1,N_d);                       % element phase
alpha_d = alpha_d*pi/180;                       % convert to radians from deg
a = 100*ones(1,N_d);                            % ???
for i=1:N_d
    Am_d(:,i) = a(i);                           % initalize amplitudes the same for all agents
    rm_d(:,i) =[(i-3)*d, 0]';                   % space elements at d distance apart with middle element resting @ origin
end


%==========================================================================
% space discretization
Num_samples = 49;                        % No. of samples from desired beam
% N_a =Num_samples+1; 
theta = linspace(0,2*pi,Num_samples);    % N samples from theta space  


%==========================================================================
% Array factor (AFd) of desired beam
AFd = GetArrayFactor(N_d, rm_d, Am_d, alpha_d, gamma, mu, nu, rho, theta); % call provided code to find AF
h = polar(theta,abs(AFd),'r'); % plot
set(h,'LineWidth',2);
box on
hold on


f = abs(AFd)';       % array factor at every sampled point
theta = theta';     
N_s =  Num_samples;  
N_a = 64;            % extend to 64 agents
rm = zeros(2,N_a);   % reinitialize positions for 64 agents
for i=1:N_a
    rm(:,i) =[(i-3)*d, 0]' ; % place 64 agents with the same spacing as before, where the 3rd is at 0
end
xm = rm(1,:)'; % split into x and y posns
ym = rm(2,:)';