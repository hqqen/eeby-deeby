close all
% xrec = rho.*cos(theta);
% yrec = rho.*sin(theta);
% scatter3(xrec,yrec,af_db(:,ind),'d', 'LineWidth', 5)
% grid on
% hold on
% set(gca, 'LineWidth', 5, 'FontSize', 35)
% xlabel('position x')
% ylabel('position y')
% zlabel('R_x power(dB)')
ind = 832;
polarplot(theta(2:end),af(2:end,ind)/max(af(:,ind)), 'LineWidth', 5)
hold on
polarplot(theta(1),af(1,ind)/max(af(:,ind)), 'o', 'LineWidth', 5)
set(gca, 'LineWidth', 5, 'FontSize', 35)
avg = mean(af(2:end,ind)/max(af(:,ind)));
dev = std(af(2:end,ind)/max(af(:,ind)));
peak = max(af(2:end,ind)/max(af(:,ind)));
disp(['avg: ',num2str(avg),' dev: ',num2str(dev), ' peak: ',num2str(peak)])