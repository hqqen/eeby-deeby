clc

%%
for m=1:N_a
    for j=1:N_s
%         p(m,j) =  gamma_ch(m,j).*exp(zeta(m,j)*1i);
        p(m,j) =  1./d(m,j)*exp(zeta(m,j)*1i);
    end
end
p = p';
pr = real(p);
pi = imag(p);
P = [pr -pi; pi pr];
rank(P)