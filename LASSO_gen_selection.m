% selecting the effective genes for estimation of each special gene and aggregate the result 
te = load('test.mat');
tr = load('train.mat');
va = load('validation.mat');
X=  [tr.x,tr.y];
allBroadData = X;
all_lasso_D=[];
all_lasso_B=[];
meanErrs=[];

for k=1:size(allBroadData,2)
    k
    temp = normc(allBroadData);
    res = temp(:,k);
    temp(:,k)=[];
    reg = temp;
    [lasso_broad,info_broad]=lasso(reg,res);
    [Ml,Il] = min(info_broad.MSE);
    result = reg*lasso_broad(:,Il)+info_broad.Intercept(Il);
    all_lasso_D = [all_lasso_D,lasso_broad(:,Il)];
    all_lasso_B = [all_lasso_B,info_broad.Intercept(Il)];
    err_lasso = result-allBroadData(:,k);   
    meanErrLasso = sqrt(sum(err_lasso.^2)/size(err_lasso,1))
    meanErrs = [meanErrs,meanErrLasso];
end
save('all_lasso_result','all_lasso_D','all_lasso_B','meanErrs');
