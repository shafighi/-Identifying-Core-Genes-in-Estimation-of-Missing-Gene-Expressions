te = load('test.mat');
tr = load('train.mat');
va = load('validation');
%X=  [tr.x,tr.y;va.x,va.y];

%%%%%%%%%%%%%%%%%%%%% parameters and selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%
landmark_size = 943;

% threshold control the maximum error which we want to accept
threshold=0.2;
genes =[1:landmark_size]; %randi(c_df-c,[1 10000]);
validation_count = size(va.x,1);
my_tr = [];

%selection of the samples which are near to the complete data
for sample=1:validation_count
    a = bsxfun(@minus,[tr.x,tr.y],[va.x(sample,:),va.y(sample,:)]);
    out = cellfun(@norm,num2cell(a,2));
    [M,I] = min(out)
    my_tr = [my_tr,[tr.x(I,:),tr.y(I,:)]]
end
X=  [my_tr;[va.x,va.y]];
save('nearX','X');
r_df=size(X,1);
c_df=size(X,2);

%%%%%%%%%%%%%%%%%%%%% mean error of matrix method %%%%%%%%%%%
  %cols = IBroad;
  df = X;
  C = df((validation_count+1):r_df,1:landmark_size);
  A = df(1:validation_count,1:landmark_size);
  %temp=df;
  df(:,1:landmark_size)=[];
  %X=[temp(:,cols),df(:,genes)];
  DX = df((validation_count+1):r_df,:);
  B = df(1:validation_count,:);
  estOfBroad = (C * pinv(A)) * B;
  errB = sum(sum(abs(estOfBroad-DX).^2));
  meanErrOfBroad = sqrt(errB / (size(DX,1)*size(DX,2)));
  badPredictionsOfBroad = (abs(estOfBroad-DX).^2)>threshold; 
  badPredictionCountBroad = sum(badPredictionsOfBroad(:) == 1);
  errorOfFirst = (estOfBroad-DX);
  save('errorOfFirst','errorOfFirst');
  all = reshape(errorOfFirst,[1,size(errorOfFirst,1)*size(errorOfFirst,2)]);
  
%%%%%%%%%%%%%%%%%%%%% Itrative algorithm %%%%%%%%%%%%%%%%%%%%%
iter = 1;
promoteBadPrediction=[];
promoteMeanError=[];
allOfMeanErr=[];
allOfBadPrediction=[];
%bestCols = con_1000;
bestCols = [1:943];%randi([1 c_df],1,c); %random initialization of columns
  cols = bestCols;
  df = X;
  C = df((validation_count+1):r_df,cols);
  A = df(1:validation_count,cols);
  df(:,cols)=[];
  DX = df((validation_count+1):r_df,:);
  B = df(1:validation_count,:);
 
  est = (C * pinv(A)) * B;
  err = sum(sum(abs(est-DX).^2));
  meanErr = sqrt(err / (size(DX,1)*size(DX,2)));
  promoteMeanError = [promoteMeanError,meanErr];
  badPredictions = (abs(est-DX).^2)>threshold; 
  badPredictionCount = sum(badPredictions(:) == 1);
  bestPredictionCount = badPredictionCount;
  promoteBadPrediction = [promoteBadPrediction,badPredictionCount];
  counter=0;
  errorOfEle =  (abs(est-DX).^2);
  errCol = sum(errorOfEle>0.02);
  w = sum(abs(pinv(A)*B),2);

for iteration = 1:10000
  df = X;
  cols = bestCols; 
  counter=counter+1;
  iteration
  [M,I] = max(errCol);
  errCol(I)=[-10000000000];
  [Lia,Locb] = ismember(df(:,I)',X','rows');
  
  
  [mw,iw] = min(w);
  w(iw)=[10000000000];
  cols(iw) = Locb;

 
  C = df((validation_count+1):r_df,cols);
  A = df(1:validation_count,cols);
  df(:,cols)=[];
  DX = df((validation_count+1):r_df,:);
  B = df(1:validation_count,:);
  est = (C * pinv(A)) * B;
  badPredictions = (abs(est-DX).^2)>threshold; 
  badPredictionCount = sum(badPredictions(:) == 1);
  allOfMeanErr = [allOfMeanErr,meanErr];
  allOfBadPrediction = [allOfBadPrediction,badPredictionCount];
    
  if badPredictionCount<bestPredictionCount
     counter=0;
     %print(iteration)
     bestPredictionCount = badPredictionCount;
     bestCols = cols;
     iter = iteration;
     bestEst = est;
     err = sum(sum(abs(est-DX).^2));
     meanErr = sqrt(err / (size(DX,1)*size(DX,2)));
     promoteMeanError = [promoteMeanError,meanErr];
     promoteBadPrediction = [promoteBadPrediction,badPredictionCount];
     errorOfEle =  (abs(est-DX).^2);
     errCol = sum(errorOfEle>0.02);
     w = sum(abs(pinv(A)*B),2);    
  end
end


%common = intersect(landmark,bestCols);
%methods = intersect(table2array(lm)',table2array(p53(bestCols,1))');
%commonGenes =p53(methods,1);
%plot(promoteBadPrediction)
%plot(promoteMeanError)
%plot(allOfBadPrediction)
%plot(allOfMeanErr)

save('first_run','allOfMeanErr','promoteBadPrediction','promoteMeanError','allOfBadPrediction','bestCols');

