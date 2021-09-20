function [left_func,right_func] = processHCPcifti(func_cifti,smoothingKernel)

	disp(['Processing ',func_cifti,' ...']);
	TMP_dir = '...'; % A path for the temp files to live in
    
	if(~isdir(TMP_dir))
		mkdir(TMP_dir)
	end
	tmp_base=[tempname(TMP_dir)];
	tmp_cifti=[tmp_base,'.dtseries.nii'];
	tmp_gifti=[tmp_base,'.func.gii'];

	% Locations of surfaces, you must have them in the right path
	left_surf=['.../Q1-Q6_RelatedParcellation210.L.midthickness_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii'];
	right_surf=['.../Q1-Q6_RelatedParcellation210.R.midthickness_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii'];
  %left_surf=['/home/ychen/kg98/Schaefer2018_400Parcels_7Networks_order.L.label.gii'];% You can change the atlas to the Schaefer atlas 
  %right_surf=['/home/ychen/kg98/Schaefer2018_400Parcels_7Networks_order.R.label.gii'];
    
	if(nargin<2)
		smoothingKernel=0;
	end

	% Do presmoothing in CIFTI space
	if(smoothingKernel>0),
		% Conversion of smoothing kernel to sigma
		sigma=smoothingKernel/sqrt(8*log(2));
		unix_command=['wb_command -cifti-smoothing ',func_cifti,' ',num2str(sigma),' ',num2str(sigma),' COLUMN ',tmp_cifti,' -left-surface ',left_surf,' -right-surface ',right_surf];
		system(unix_command);
	else
		% No smoothing option
		tmp_cifti=func_cifti;
	end
	
	% Now take the cifti and convert it to a gifti
	unix_command=['wb_command -cifti-convert -to-gifti-ext ',tmp_cifti,' ',tmp_gifti];
	system(unix_command);
	% disp(tmp_gifti);
	% Then after this has been done, now load in the cifti into matlab
 
    gifti_loaded=gifti(tmp_gifti);
	[left_func,right_func] = fromGrayOrdinatesToSurface(gifti_loaded);
	% Next step is to delete all the temporary files

	% Now do the transformation to look at precent signal change
	left_func=(left_func-repmat(mean(left_func,2),[1 size(left_func,2)]))./repmat(mean(left_func,2),[1 size(left_func,2)]);
	right_func=(right_func-repmat(mean(right_func,2),[1 size(right_func,2)]))./repmat(mean(right_func,2),[1 size(right_func,2)]);

	% Deal with nans coming from the division:
	m_l=mean(left_func,2);left_func(isnan(m_l),:)=0;
	m_r=mean(right_func,2);right_func(isnan(m_r),:)=0;

	unix_command=['rm -rf ',tmp_base,'.*'];
	system(unix_command);
end
