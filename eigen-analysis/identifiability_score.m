function[maxIDs,maxScale,IDScore] = identifiability_score(mat_t1,mat_t2,max_iter)
%% This function calculates the identifiability score 
%% between the eigenvalue matrixes of time 1 (mat_t1) and time 2 (mat_t2).
%% Both mat_t1 and mat_t2 are subjects*eigenvalues(or can be whatever measurements: 
%% e.g., each column for thickness of each region);
%% mat_t1 and mat_t2 should have the same size and have at least 2 columns
%% Without specifying max_iter, the identifiability score(IDScore) is a scaler, 
%% which is calculated from the whole matrices (mat_t1, mat_t2)
%% Setting the max_iter to iteratively running the identifiability score from the first 
%% 2 columns of the matrices to the first max_iter columns, e.g., 1000
%% maxIDs is the maxinum identifiability score across the first max_iter columns, if max_iter is specified.
%% maxEV is the column (or eigenvalue) index corresponding to the maxIDs   

no_sub=size(mat_t1,1);
idx = find(~eye(size(zeros(no_sub))));

if nargin == 2 % running the identifiability score for the whole matrices
    IDScore=zeros(1);
    
    [r,~] = corr(mat_t1', mat_t2');
    x1m = mean(diag(r)); %x1m is the mean of within-subject correlation
    
    x2m = mean(r(idx)); %x2m is the mean of between-subject correlation
    s2 = std(r(idx)); %x2m is the std of between-subject correlation
    
    IDScore = (x1m-x2m)/s2;
    maxIDs = []; maxScale = [];
end

if nargin == 3 % running the identifiability score iteratively
    IDScore = zeros(1,max_iter);
    
    for scale = 1:max_iter;
        [r,~] = corr(mat_t1(:,1:scale)', mat_t2(:,1:scale)');
        x1m = mean(diag(r)); %x1m is the mean of within-subject correlation
        
        x2m = mean(r(idx)); %x2m is the mean of between-subject correlation
        s2 = std(r(idx)); %x2m is the std of between-subject correlation
        
        IDScore(scale) = (x1m-x2m)/s2;
    end
    
    maxScale = find(IDScore == max(IDScore));
    maxIDs = max(IDScore);
    %figure; plot(IDScore);
    %xlabel('Eigenvalue Index')
    %ylabel('Identifiability Score')
end
end
