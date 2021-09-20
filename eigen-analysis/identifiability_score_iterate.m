function[maxIDs,maxEV,IDScore_it] =identifiability_score_iterate(EV_t1,EV_t2,EV_no)
%% This function calculates the identifiability score iteratively across a range of eigenvalues
%% between the eigenvalue matrixes of time 1 (EV_t1) and time 2 (EV_t2).
%% Both EV_t1 and EV_t2 are subjects*eigenvalues; EV_no is the maxinum number of eigenvalues, e.g., 1000
%% The IDScore_it is a vector of the identifiability score; maxIDs is the maxinum identifiability score.
%% maxEV is the eigenvalue index corresponding to the maxIDs   
%% If you'd like to calculate the identifiability score based on single iteration 
%% (e.g., just for the first 144 eigenvalues), plese check "identifiability_score_index.m" 

no_sub=size(EV_t1,1);
IDScore_it=zeros(1,EV_no);
idx = find(~eye(size(zeros(no_sub))));

for EV=1:EV_no;
    [r,p]=corr(EV_t1(:,1:EV)', EV_t2(:,1:EV)');
    x1m=mean(diag(r)); %x1m is the mean of within-subject correlation
  
    x2m=mean(r(idx)); %x2m is the mean of between-subject correlation
    s2=std(r(idx)); %x2m is the std of between-subject correlation
    
    IDScore_it(EV)=(x1m-x2m)/s2;
end

maxEV=find(IDScore_it==max(IDScore_it));
maxIDs=max(IDScore_it);
%figure;plot(IDScore);
%xlabel('Eigenvalue Index')
%ylabel('Identifiability Score')
