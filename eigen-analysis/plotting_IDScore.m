figure('color','white');
plot(IDScore_oasis_surf,'lineWidth',2); %can be changed to other dataset, e.g., adni or hcp
%xlabel('Eigenvalue Index')
%ylabel('Identification Score')
set(gca,'fontSize',24,'fontname','Arial','box','on')
hold on
plot(IDScore_oasis_surf_L,'lineWidth',2);
plot(IDScore_oasis_surf_R,'lineWidth',2);
plot(IDScore_oasis_surf_LR_con,'lineWidth',2);
plot(IDScore_oasis_surf_LtoR_T1,'b--','lineWidth',2);
plot(IDScore_oasis_surf_LtoR_T2,'r--','lineWidth',2);
yline(mean(maxIDs_oasis_perm_new_all),'color','k','lineWidth',2);
ymin=mean(maxIDs_oasis_perm_new_all)-(2*std(maxIDs_oasis_perm_new_all)); 
ymax=mean(maxIDs_oasis_perm_new_all)+(2*std(maxIDs_oasis_perm_new_all))
y1=yline(ymax,'color','[0.95,0.95,0.95]','lineWidth',0.001);
y2=yline(ymin,'color','[0.95,0.95,0.95]','lineWidth',0.001);
patch([0 1000 1000 0], [ymin ymin, ymax ymax], [0.55 0.55 0.55],'EdgeColor','none');
alpha(0.4)% alpha set patch transparency rate
y1.Color(4) = 0.1;
y2.Color(4) = 0.1;
xticks([0:200:1000])
ylim([-1 7])
yticks([-1,0,1,3,5,7])
%exportgraphics(gcf,'IDS_oasis_all.jpg','Resolution',900)
