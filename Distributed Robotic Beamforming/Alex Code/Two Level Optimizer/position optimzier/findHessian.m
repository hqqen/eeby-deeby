syms wi fi am gammami dmi alpham k xm ym thti rhoi real
d = norm([xm; ym] - rhoi*[cos(thti); sin(thti)]);
L = (wi/2)*norm(fi-norm((am*gammami/dmi)*exp(1j*(alpham + k*(xm*cos(thti) + ym*sin(thti) + dmi)))));

lGrad = [diff(L,xm); diff(L,ym)];
lHessian = hessian(L,[xm,ym]);