%========================================================================= 
clc
clear all
close all
rng default

%==========================================================================
% Parameters of desired beam
f = 40e6;                                     % Frequency in Hz                                                                        
N_d =5;                                       % Number of elements of array  
lambda = 3e8/f;                               % wave length
d = lambda/4;                                 % Inter-element spacing
k = (2*pi)./lambda;
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
alpha_d = alpha_d*pi/180;
a = 100*ones(1,N_d);
for i=1:N_d
    Am_d(:,i) = a(i);
    rm_d(:,i) =[(i-3)*d, 0]' ;
end


%==========================================================================
% space discretization
Num_samples = 49;                        % No. of samples from desired beam
% N_a =Num_samples+1; 
theta =linspace(0,2*pi,Num_samples);     % N samples from theta space  


%==========================================================================
% Array factor (AFd) of desired beam
AFd = GetArrayFactor(N_d, rm_d, Am_d, alpha_d, gamma, mu, nu, rho, theta);
h=polar(theta,abs(AFd),'r');
set(h,'LineWidth',2);
box on
hold on


f = abs(AFd)';
theta = theta';
N_s =  Num_samples;
N_a = 64;
rm = zeros(2,N_a);  
for i=1:N_a
    rm(:,i) =[(i-3)*d, 0]' ;
end
xm = rm(1,:)';
ym = rm(2,:)';