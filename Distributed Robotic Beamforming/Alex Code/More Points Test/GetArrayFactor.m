%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Name: GetArrayFactor.m
%% Description: This function computes  array factor
%% Inputs: N (number of elements)  rm (element's location)  Am
%%         (amplitude)  alpha (phase)  gamma, mu, nu (ch. parameters)
%% Output:   AF
%% SubFunctions:    None
%% Last modified: 11/13/2019 by Anjaly Parayil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AF = GetArrayFactor(N, rm, Am, alpha, gamma, mu, nu, rho, theta)


%wave parameters
f = 40e6;                           % Frequency in Hz
lambda = 3e8/f;
d = .25*lambda;%/4;                       % Inter-element spacing

%Derived parameters
fu = f*ones(1,N);                           %frequency vector with errors
l = 3e8./fu;                                %wavelength in meters with errors
k = (2*pi)./l;                              %wavenumber with frequency uncertainity


%% Compute array factor (AF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AF = zeros(length(rho),length(theta));
r_hat = [cos(theta); sin(theta)];       % farfiled directions
ch=[1+2i;2-2i;3-i;1;1];
for ind = 1:length(rho)
    af = zeros(1,length(theta));


    for kk = 1:N
        mm = kk - 1;
        Im = Am(1,kk)*exp(1i*alpha(1,kk));    % element excitation
        dm_mat = (kron(ones(1,length(theta)),rm(:,kk)) - rho(ind)*r_hat);
        dm = sqrt(sum(dm_mat.^2, 1)) + 1;     % distance between element and receiver
        %        ch(kk,:) = gamma./((dm).^(mu/2));     % channel


        % with channel model
        psi1 =  k(1,kk)*rm(:,kk)'*r_hat ;%+ nu*k(1,kk)*dm;
        af = af + Im*exp(1i*psi1);%.*ch(kk,:);

    end

    AF(ind,:) = af;
end


%
% for ind = 1:length(rho)
%     af1 = zeros(N,length(theta));
%
%
% for kk = 1:N
%         mm = kk - 1;
%         Im = Am(1,kk)*exp(1i*alpha(1,kk));    % element excitation
%         dm_mat = ( kron(ones(1,length(theta)),rm(:,kk)) - rho(ind)*r_hat );
%         dm = sqrt(sum(dm_mat.^2, 1)) + 1;     % distance between element and receiver
%         ch(kk,:) = gamma./((dm).^(mu/2));     % channel
%
%
%         % with channel model
%         psi1 =  k(1,kk)*rm(:,kk)'*r_hat ;%+ nu*k(1,kk)*dm;
%         af1(kk,:) =  N*Im*exp(1i*psi1);%.*ch(kk,:);
% end
% end
% keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

