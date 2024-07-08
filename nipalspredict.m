function mdl=nipalspredict(X,mdl,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                           nipals predict
%                     nipals algorithm reported in:
%               https://doi.org/10.1002/cem.1180020306
%                   Paper author: A. Hoskuldsson
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% this function predict variables based on the parameters extracte by a
% classical PLS model
% - X matrix
% X matrix to analyze
% - mdl
% clasical PLS model built with the sister function
%-----------------------------------------------------------------<\header>
defaultisScaled=true;
defaultY=0;
checkInput=inputParser;
addRequired(checkInput,'X')
addRequired(checkInput,'mdl')
addParameter(checkInput,'isScaled',defaultisScaled,@islogical)
addParameter(checkInput,'Y',defaultY,@isnumeric)
parse(checkInput,X,mdl,varargin{:})
%init
X=checkInput.Results.X;
mdl=checkInput.Results.mdl;
Y=checkInput.Results.Y;
isScaled=checkInput.Results.isScaled;
%--------------------------------------------------------------------------
warning('off')
% Autoscaling (X and Y)
if isScaled
    
elseif ~isCentered
        X=normalize(X,'center',mdl.xCenter,'scale',1);
else
    X=normalize(X,'center',mdl.xCenter,'scale',mdl.xScale);
end
Xorig=X;
nComp=size(mdl.T,2);
T_predict=nan(size(X,1),nComp);
U_predict=nan(size(X,1),nComp);
for i=1:nComp
    T_predict(:,i)=X*mdl.W(:,i)/(mdl.P(:,i)'*mdl.W(:,i));
    X = X - T_predict(:,i)*mdl.P(:,i)';
    U_predict(:,i)=mdl.B(:,i)*T_predict(:,i);
end
Yhat                =U_predict*mdl.Q';
Ypredict            =Yhat;
mm                  =size(mdl.T,1);
kk                  =size(mdl.P,1);
Eig                 =diag(diag((mdl.T'*mdl.T)/(mm-1)));
Ht2predict          =diag(T_predict/(Eig)*T_predict');
SPEpredict          =diag(Xorig*(eye(kk)-orth(mdl.P)*orth(mdl.P)')*Xorig');
predict.Y           =Ypredict;
predict.X           =X;
predict.T           =T_predict;
predict.T_s         =[];
predict.U           =U_predict;
predict.W           =mdl.W;
predict.R2Y         =nan;
predict.RMSE        =nan;
predict.MAE         =nan;
predict.Hotelling   =Ht2predict;
predict.SPE         =SPEpredict;
predict.Ytrue       =nan(size(Ypredict));
if all(Y==0)
else
    nanMask=~any(ismissing(U_predict),2);
    TSS=    mean((Y(nanMask,:)).^2,'all','omitnan');%-mean(Y(nanMask,:))
    PRESS=  mean((Y-U_predict.*mdl.Q).^2,1,'omitnan');
    WoldRCriterion=PRESS(2:end)./PRESS(1:end-1);
    Q2Y=    1-(PRESS/TSS);
    R2Y=    nan(1,size(Y,2));
    RMSE=   nan(1,size(Y,2));
    MAE=    nan(1,size(Y,2));
    for i=1:size(Y,2)
        ESS=        mean((Y(:,i)-Ypredict(:,i)).^2,1,'omitnan');
        R2Y(i)=     1-(ESS/TSS);
%         R2Y(i)=     corr(Y(:,i),Ypredict(:,i),'Rows','complete');
        RMSE(i)=    sqrt(mean((Y(:,i)-Ypredict(:,i)).^2,'omitnan'));
        MAE(i)=     mean(abs(Y(:,i)-Ypredict(:,i)),'omitnan');
    end
    predict.Q2Y=Q2Y;
    predict.R2Y=R2Y;%.^2;
    predict.Ytrue=Y;
    predict.RMSE=RMSE;
    predict.MAE=MAE;
    predict.WoldRCriterion=WoldRCriterion;
end
mdl.predict=predict;
