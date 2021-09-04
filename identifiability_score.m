function[maxIDs,maxEV,IDScore] =identifiability_score(EV_T1,EV_T2,EV_no)
%% This function calculates the identifiability score across a range of eigenvalues
%% between the eigenvalue matrixes of time 1 (EV_T1) and time 2 (EV_T2).
%% Both EV_T1 and EV_T2 are subjects*eigenvalues; EV_no is the number of eigenvalues, e.g., 1000
%% The IDScore is a vector of the identifiability score; maxIDs is the maxinum identifiability score.
%% maxEV is the eigenvalue index corresponding to the maxIDs   

no_sub=size(EV_T1,1);
IDScore=zeros(1,EV_no);
idx = find(~eye(size(zeros(no_sub))));

for i=1:EV_no;
    [r,p]=corr(EV_T1(:,1:i)', EV_T2(:,1:i)');
    x1m=mean(diag(r)); %x1m is the mean of within-subject correlation
  
    x2m=mean(r(idx)); %x2m is the mean of between-subject correlation
    s2=std(r(idx)); %x2m is the std of between-subject correlation
    
    IDScore(i)=(x1m-x2m)/s2;
end
%figure;plot(IDScore);
maxEV=find(IDScore==max(IDScore));
maxIDs=max(IDScore);
