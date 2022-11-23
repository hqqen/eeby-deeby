syms q1 q2 q3 q4 q5 q6 d1 d2 d3 d4 d5 d6 a1 a2 a3 a4 a5 a6 alp1 alp2 alp3 alp4 alp5 alp6

DH = [...
      q1, d1, a1, alp1;
      q2, d2, a2, alp2;
      q3, d3, a3, alp3;
      q4, d4, a4, alp4;
      q5, d5, a5, alp5;
      q6, d6, a6, alp6];


T = []
for i = 1:max(size(DH))

    tht = DH(i,1); d = DH(i,2); a = DH(i,3); alp = DH(i,4);



    H(i,:,:) = hessian(T,q);

end
