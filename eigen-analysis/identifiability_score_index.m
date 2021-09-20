function[IDScore]=identifiability_score_index(mat_t1,mat_t2);
%% This function calculate the identifiability score of the two matrices: mat_t1 and mat_t2
%% The row of mat_t1 and mat_t2 is subject; both matrices should have at least two columns 
no_sub=size(mat_t1,1);
IDScore=zeros(1);
idx = find(~eye(size(zeros(no_sub))));

[r,p]=corr(mat_t1', mat_t2');
x1m=mean(diag(r)); %x1m is the mean of within-subject correlation
  
x2m=mean(r(idx)); %x2m is the mean of between-subject correlation
s2=std(r(idx)); %x2m is the std of between-subject correlation
    
IDScore(i)=(x1m-x2m)/s2;
end
