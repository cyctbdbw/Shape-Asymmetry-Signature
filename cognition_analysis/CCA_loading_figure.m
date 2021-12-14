[r_sas,~]=corr(U(:,1),sas_mean); % U is the SAS canonical variate; 
% sas_mean is the mean of SAS across each eigen-group
figure('color','white');
b=bar(r_sas,'b');
b.FaceColor = 'flat';
un_sig_idx=find(padj_sas>=0.05); 
% padj_sas is FDR corrected p value of the CCA loading, 
% which is calculated by bootstrapping using the CCA_sig_test.R

for id=1:length(un_sig_idx)
    b.CData(un_sig_idx(id),:) = [0.85 0.85 0.85];
end

hold on

er = errorbar([1:length(r_sas)],r_sas,se_sas,'LineWidth',1,'CapSize',20); 
%se_sas is from the bootstrapping

er.LineStyle = 'none';
er.Color = [0 0 0];
set(gca,'fontSize',20,'box','on')
a=padj_sas';
xlabel('Eigen-group')
ylabel('Loading')
yticks([-1:0.2:1])
ylim([-1 1])
exportgraphics(gcf,'SAS_CCA_loading.jpg','Resolution',900)
