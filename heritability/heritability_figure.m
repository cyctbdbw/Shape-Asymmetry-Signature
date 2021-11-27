TL=readtable('heritabilityACTE_white_L_all_resultsNA.txt'); %
TL=table2array(TL);
pL=TL(:,8);

[padj_L,~] = multicmp (pL(2:end),'fdr');
sig_L_ind=find(padj_L<0.05);
sig_L_ind_new=sig_L_ind+1;

figure('color','white');
b=bar([0;TL(2:end,1)],'FaceColor','[0.85 0.85 0.85]','EdgeColor','none'); %We only show the 2-144 eigenvalues
b.FaceColor = 'flat';
for i=1:length(sig_L_ind_new)
       b.CData(sig_L_ind_new(i),:) = [0.6350, 0.0780, 0.1840];
end
%ylabel('Heritability')
%xlabel('Eigenvalue Index')
set(gca,'fontSize',26,'box','on','fontname','Arial')
xlim([0 145]);
ylim([0 1])
%title 'Left Hemisphere'
xtl ={{'0';''},{'';'1'},{'';'3'},{'20';''},{'';'4'},{'';'5'},{'';'6'},{'50';''},{'';'7'},{'';'8'},{'';'9'},{'100', ''},{'';'10'},{'';'11'},{'144';''}};
h = my_xticklabels(gca,[0,3,13,20,21,31,43,50,57,73,91,100,111,133,144],xtl, 'FontSize',26,'fontname','Arial');
yticks([0.2 0.4 0.6 0.8 1]);
xline(144,'--','color','k','lineWidth',1);
xline(121,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(100,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(81,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(64,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(49,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(36,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(25,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(16,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(9,'--','color',[0.7,0.7,0.7],'lineWidth',1);
xline(4,'--','color',[0.7,0.7,0.7],'lineWidth',1);
%exportgraphics(gcf,'Heritibility_surf_L.jpg')
