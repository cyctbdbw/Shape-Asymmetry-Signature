parcellation='Schaefer' %You need to choose Schaefer or MMP
TimePoint='T1' %You need to choose T1 or T2
Schaefer_N=100; %You need to choose the N if the parcellation is 'Schaefer'

switch TimePoint
    case 'T1'
        cd /home/ychen/kg98/ychen/HCP_T1_fMRI;
        baseFolder='.../HCP_T1_fMRI';
        % You need to chage the path where the HCP_T1_fMRI live
    case 'T2'
        cd /home/ychen/kg98_scratch/Yuchi/HCP_T2_fMRI;
        baseFolder='.../HCP_T2_fMRI';
        % You need to chage the path where the HCP_T2_fMRI live
end

load('subj.mat'); %This is the subject list

smoothingKernel=0; % no smoothing; you can change the smoothing level.

switch parcellation
    case 'MMP'
        labelGifti_left=gifti('Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii');
        labelGifti_right=gifti('Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.label.gii');
        parcel_flag=360 %for the HCP
        uT=find(triu(ones(360,360),1));
    case 'Schaefer'
        labelGifti_left=gifti(['Schaefer2018_',num2str(Schaefer_N),'Parcels_7Networks_order.L.label.gii']);
        labelGifti_right=gifti(['Schaefer2018_',num2str(Schaefer_N),'Parcels_7Networks_order.R.label.gii']);
        uT=find(triu(ones(Schaefer_N,Schaefer_N),1));
        % % Below is a fix to replace the gifti with a structure to make the same parcellation code work.
        cdata=labelGifti_left.cdata;labels=labelGifti_left.labels;
        labelGifti_left=struct;labelGifti_left.labels=labels;labelGifti_left.cdata=cdata;
        cdata=labelGifti_right.cdata;labels=labelGifti_right.labels;
        labelGifti_right=struct;labelGifti_right.cdata=cdata;labelGifti_right.labels=labels;
        
        % % Seperate left and right
        labelGifti_left.labels.key = labelGifti_left.labels.key([1:(Schaefer_N/2+1)]);
        labelGifti_left.labels.name = labelGifti_left.labels.name([1:(Schaefer_N/2+1)]);
        labelGifti_right.labels.key = labelGifti_right.labels.key([1,(Schaefer_N/2+2):(Schaefer_N+1)]);
        labelGifti_right.labels.name = labelGifti_right.labels.name([1,(Schaefer_N/2+2):(Schaefer_N+1)]);
end

parcel_flag=360;

for subject=1:45; %45 is the total number of subjects
    try
        subject_name=num2str(subj{subject});
        
        subject_folder_LR=[baseFolder,'/',subject_name,'/MNINonLinear/Results/rfMRI_REST1_LR'];
        subject_folder_RL=[baseFolder,'/',subject_name,'/MNINonLinear/Results/rfMRI_REST1_RL'];
        
        switch parcellation
            case 'MMP'
                [left_func_1LR,right_func_1LR] = processHCPcifti([subject_folder_LR,'/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'],smoothingKernel);
                [left_func_1RL,right_func_1RL] = processHCPcifti([subject_folder_RL,'/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'],smoothingKernel);
            case 'Schaefer'
                [left_func_1LR,right_func_1LR] = processHCPcifti_Schaefer([subject_folder_LR,'/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'],smoothingKernel,Schaefer_N);
                [left_func_1RL,right_func_1RL] = processHCPcifti_Schaefer([subject_folder_RL,'/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'],smoothingKernel,Schaefer_N);
        end
        
        % For the FCMs we need to save them
        time1_left=[left_func_1LR,left_func_1RL];
        time1_right=[right_func_1LR,right_func_1RL];
        clear left_func_1RL left_func_1LR right_func_1LR right_func_1RL
        
        switch parcellation
            case 'MMP'
                time1_original(:,:)=[parcelTimeSeries_Schaefer(time1_left,labelGifti_left,360);parcelTimeSeries_Schaefer(time1_right,labelGifti_right,360)];
            case 'Schaefer'
                time1_original(:,:)=[parcelTimeSeries_Schaefer(time1_left,labelGifti_left,Schaefer_N);parcelTimeSeries_Schaefer(time1_right,labelGifti_right,Schaefer_N)];
        end
        
        [r,p]=corr(time1_original');
        
        FC_all=triu(r,1);
        FC=FC_all(find(FC_all));,
        
        switch parcellation
            case 'MMP'
                switch TimePoint
                    case 'T1'
                        save([subject_name,'_FC_a_MMP'],'FC');
                    case 'T2'
                        save([subject_name,'_FC_b_MMP'],'FC');
                end
            case 'Schaefer'
                switch TimePoint
                    case 'T1'
                        save([subject_name,'_FC_a_',num2str(Schaefer_N)],'FC');
                    case 'T2'
                        save([subject_name,'_FC_b_',num2str(Schaefer_N)],'FC');
                end
        end
        
    catch
        disp(['no ',subject_name])
    end
end
