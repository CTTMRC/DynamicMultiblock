function mdl=DiPLSpredict(X,mdl,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                           DiPLS predict
%                     DiPLS algorithm reported in:
%           http://dx.doi.org/10.1016/j.ifacol.2015.08.167
%                  Paper authors: Yining Dong and S. Joe Qin
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% this function predict variables based on the parameters extracte by a
% DiPLS model
% - X matrix
% X matrix to analyze
% - mdl
% DiPLS model built with the sister function
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
% Autoscaling (X and Y)
if isScaled
    
elseif ~isCentered
    X=normalize(X,'center',mdl.xCenter,'scale',1);
else
    X=normalize(X,'center',mdl.xCenter,'scale',mdl.xScale);
end
nComp=size(mdl.T,2);
s=size(mdl.BetaX,1)-1;
T_predict=nan(size(X,1),nComp);
U_predict=nan(size(X,1)-s,nComp);
for i=1:nComp
    T_predict(:,i)=X*mdl.W(:,i);
    X = X - T_predict(:,i)*mdl.P(:,i)';
    T_s_predict=nan(size(X,1)-s,s+1);
    T_s_predict(:,1,i)=T_predict(s+1:end,i);
    for j=1:s
        T_s_predict(:,j+1,i)=T_predict(s-j+1:end-j,i);
    end
    U_predict(:,i)=T_s_predict(:,:,i)*mdl.BetaXhat(:,i);
end
Yhat                =U_predict*mdl.Q';
Ypredict            =[nan(s,size(Yhat,2));Yhat];
mm                  =size(mdl.T,1);
kk                  =size(mdl.P,1);
Eig                 =diag(diag((mdl.T'*mdl.T)/(mm-1)));
Ht2predict          =diag(X*orth(mdl.P)/(Eig)*orth(mdl.P)'*X');
SPEpredict          =diag(X*(eye(kk)-orth(mdl.P)*orth(mdl.P)')*X');
predict.Y           =Ypredict;
predict.X           =X;
predict.T           =T_predict;
predict.T_s         =T_s_predict;
predict.U           =U_predict;
predict.W           =mdl.W;
predict.R2Y         =nan;
predict.RMSE        =nan;
predict.Hotelling   =Ht2predict;
predict.SPE         =SPEpredict;
if all(Y==0)
else
    
    TSS=    mean((Y(s+1:end,:)).^2,1,'omitnan');%-mean(Y(nanMask,:))
    PRESS=  mean((Y(s+1:end,:)-U_predict.*mdl.Q).^2,1,'omitnan');
    WoldRCriterion=PRESS(2:end)./PRESS(1:end-1);
    Q2Y=    1-(PRESS/TSS);
    R2Y=    nan(1,size(Y,2));
    RMSE=   nan(1,size(Y,2));
    MAE=    nan(1,size(Y,2));
    for i=1:size(Y,2)
        ESS=        mean((Y(:,i)-Ypredict(:,i)).^2,1,'omitnan');
        R2Y(i)=     1-(ESS/TSS);
%         R2Y(i)=     corr(Y(:,i),Ypredict(:,i),'Rows','complete');
        RMSE(i)=    sqrt(mean((Y(s+1:end,i)-Ypredict(s+1:end,i)).^2,'omitnan'));
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
