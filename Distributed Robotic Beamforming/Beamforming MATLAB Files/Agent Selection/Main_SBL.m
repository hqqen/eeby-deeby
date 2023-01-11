%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Name: Main_SBL.m
%% Description: This program implements sparse bayesian learning to select 
%%              a subset of agents from a given set in order to resynthesize 
%%              a given beampattern
%%              
%% Inputs: set of N samples from desired beam (correspond to theta values)
%%         AFd represents array factor of desired beam
%% Output: Selected subset of agents, AF: Array factor of actual beam
%% Subfiles: GetArrayFactor.m: returns desired array factor
%%%%%%% MSBL.m: returns real and imaginary part of excitation (w) for each agent
%% Last modified: by Anjaly Parayil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [] = Main_SBL()

%========================================================================= 
clc
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
N_a =Num_samples+1; 
theta =linspace(0,2*pi,Num_samples);     % N samples from theta space 


% Array factor (AFd) of desired beam
AFd = GetArrayFactor(N_d, rm_d, Am_d, alpha_d, gamma, mu, nu, rho, theta);
h=polar(theta,real(AFd),'r');
set(h,'LineWidth',2);
box on


% True parameters (i.e., initial conditions)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 100*ones(1,N_a);

Am = abs(a + 10*randn(1,N_a));                        % element amplitude
rm = zeros(2,N_a);                                    % element locations
alpha =45*ones(1,N_a) + 3*randn(1,N_a);                % element phase
alpha = alpha*pi/180;
beta = 45*pi/180;                                   % rotate the linear array by beta
Rmat = [cos(beta) -sin(beta);sin(beta) cos(beta)];  % Rotation matrix
for i=1:N_a %:-1:1
    pos=i;N_a-i-1;
    rm(:,i) =+1*randn(2,1);2*[0.013517; 0.051525];%
   
end

                                   
% The below part is to separate the real and imaginary part of the phi matrix
%==========================================================================
H=zeros(length(AFd),N_a);
phi_1=zeros(Num_samples,N_a);
phi_2=zeros(Num_samples,N_a);
% projection matrix
for l=1:length(AFd)
for agnt=1:N_a
    H(l,agnt)=exp(1i*k*(rm(1,(agnt))*cos(theta(l))+rm(2,(agnt))*sin(theta(l))));
    phi_1(l,agnt) = cos(k*(rm(1,agnt)*cos(theta(l))+rm(2,agnt)*sin(theta(l))));
    phi_2(l,agnt) = sin(k*(rm(1,agnt)*cos(theta(l))+rm(2,agnt)*sin(theta(l))));
end
end
phi=[phi_1,-phi_2;phi_2,phi_1];
Phi= phi;

% Computes real and imaginary part of excitation (X: vector of real and imaginary part of w ), lambda_SBL:  helps to select subset of agents
%==========================================================================
Y=[real(AFd');imag(AFd')];
lambda_SBL=70;
Learn_Lambda=0;
n_w=2*N_a;
for cc=1:1
[X,gamma_ind,gamma_est,count] = MSBL(Phi,Y,lambda_SBL,Learn_Lambda,'prune_gamma',1e-4,'max_iters',500,'epsilon',1e-8,'print',0);
[ind]=find(X);
ind_len=length(ind);
clear phi_H_up
phi_H_up=zeros(2*Num_samples,ind_len);
X_new=zeros(ind_len,1);
for i=1:ind_len
    phi_H_up(:,i)=Phi(:,ind(i));
    X_new(i)=X(ind(i));
end
% phi_H_up= onlinedict(Y,phi_H_up,X_new);
clear Phi
Phi=zeros(2*Num_samples,ind_len);
Phi=phi_H_up;
 
%  Below part finds array factor based on obtained w values
%==========================================================================
AFF_n2=Phi*X;
AF_up=AFF_n2(1:Num_samples)+1i*AFF_n2(Num_samples+1:2*Num_samples);
phi_H_up1=Phi(1:Num_samples,1:.5*ind_len)+1i*Phi(Num_samples+1:2*Num_samples,1:.5*ind_len);

AFF_n1=phi_H_up1*(X(1:.5*ind_len)+1i*X(.5*ind_len+1:ind_len));
norm(abs(AFF_n1)-abs(AFd)')
polar(theta,real(AFF_n1)','g')

KK=X(1:.5*ind_len)+1i*X(.5*ind_len+1:ind_len);
for agnt=1:.5*ind_len
if (abs(KK(agnt))<1)
    KK(agnt)=0;
end 
end
nnz(KK)  % output number of pruned agents
end
