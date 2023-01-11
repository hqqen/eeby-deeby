a = x(1:N_a,1);
alpha = x(N_a+1:2*N_a,1);
for m=1:N_a
    for i=1:N_s
        u_ch(m,i) = gamma_amp(m,i)*cos(alpha(m,:)+gamma_ph(m,i));
        v_ch(m,i) = gamma_amp(m,i)*sin(alpha(m,:)+gamma_ph(m,i));
    end
end
for i=1:N_s
    den_ch = 0;
    for m=1:N_a
        den_ch = den_ch + a(m,:)*(u_ch(m,i)+1i*v_ch(m,i));
    end
    den_ch = abs(den_ch);
    af_ch(i,1) = den_ch;
end
af_db(:,1) = 20*log10(af_ch(:,1)/max(af_ch(:,1)));
disp([' Rx1: ',num2str(af_db(1,1)),' Rx2: ',num2str(af_db(2,1))])