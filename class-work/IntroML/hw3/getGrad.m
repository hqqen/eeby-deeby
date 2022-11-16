syms vu vm bu bm mu y ratings lambda m u

loss = .5*(vu.'*vm + bu + bm + mu - y).^2;
% loss = symsum(loss,m,0,ratings);
% loss = symsum(loss,u,0,ratings);

penalty = (lambda/2)*(symsum(norm(vu,2)^2 + bu^2,u) + symsum(norm(vm,2)^2 + bm^2,m));

cost = loss + penalty;