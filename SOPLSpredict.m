function mdl=SOPLSpredict(X,mdl,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                           SOPLS predict
%                     SOPLS algorithm reported in:
%           https://doi.org/10.1016/j.chemolab.2007.04.002
%        Paper authors: K.Jørgensen, B.H.Mevik and T.Næs(T.Naes)
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% this function predict variables based on the parameters extracte by a
% SOPLS model
% - X cell
% cell with different X matrices to analyze
% - mdl
% SOPLS model built with the sister algorithm
%-----------------------------------------------------------------<\header>
defaultisScaled=true;
defaultY=0;
checkInput=inputParser;
addRequired(checkInput, 'X')
addRequired(checkInput, 'mdl')
addParameter(checkInput,'isScaled',defaultisScaled,@islogical)
addParameter(checkInput,'isCentered',defaultisScaled,@islogical)
addParameter(checkInput,'Y',defaultY,@isnumeric)
parse(checkInput,X,mdl,varargin{:})
%init
X=checkInput.Results.X;
mdl=checkInput.Results.mdl;
Y=checkInput.Results.Y;
isScaled=checkInput.Results.isScaled;
isCentered=checkInput.Results.isCentered;
nBlocks=length(mdl);
mType=mdl.type;
%--------------------------------------------------------------------------
% Autoscaling (X and Y)
if isScaled

elseif ~isCentered
    for i=1:nBlocks
        X{i}=normalize(X{i},'center',mdl(i).xCenter,'scale',1);
    end
    Y=normalize(Y,'center',mdl(1).yCenter,'scale',1);
else
    for i=1:nBlocks
        X{i}=normalize(X{i},'center',mdl(i).xCenter,'scale',mdl(i).xScale);
    end
    Y=normalize(Y,'center',mdl(1).yCenter,'scale',mdl(1).yScale);
end

predictY=zeros(size(X{1},1),size(mdl(1).Y_i,2),nBlocks);

switch mType
    case 'nipals'
        U_predict=nan(size(X{1},1),size(mdl(1).U,2),nBlocks);
        C=zeros(size(mdl(1).C,1),size(mdl(1).C,2),nBlocks);
        for i=1:nBlocks
            XCurr=X{i};
            if i>1
                previousBlocks=zeros(1,i-1);
                for j=1:i-1
                    previousBlocks(j)=size(mdl(j).predict.T,2);
                end

                Torth=nan(size(XCurr,1),sum(previousBlocks));
                spanStart=[0,cumsum(previousBlocks)]'+1;
                spanEnd=cumsum(previousBlocks);
                for j=1:i-1
                    nanMask=~any(ismissing(mdl(j).predict.T),2);
                    Torth(nanMask,spanStart(j):spanEnd(j))=mdl(j).predict.T(nanMask,:);
                end
                nanMask=~any(ismissing(Torth),2);
                XCurr(nanMask,:)=XCurr(nanMask,:)-Torth(nanMask,:)...
                    *pinv(Torth(nanMask,:)'*Torth(nanMask,:))*Torth(nanMask,:)'*XCurr(nanMask,:);
            end
            mdl(i)=nipalspredict(XCurr,mdl(i),'Y',Y);
            predictY(:,:,i)=mdl(i).predict.Y;
            U_predict(1:size(mdl(i).predict.U,1),1:size(mdl(i).predict.U,2),i)=mdl(i).predict.U;
            C(1:size(mdl(i).C,1),1:size(mdl(i).C,2),i)=mdl(i).C;
        end
        Ypredict=sum(predictY,3);
    case 'DiPLS'
        U_predict=nan(size(X{1},1),size(mdl(1).U,2),nBlocks);
        C=zeros(size(mdl(1).B,1),size(mdl(1).B,2),nBlocks);
        for i=1:nBlocks
            XCurr=X{i};
            %             nanMask=~any(ismissing(mdl(i).T),2);
            if i>1
                previousBlocks=zeros(1,i-1);
                for j=1:i-1
                    previousBlocks(j)=size(mdl(j).predict.T,2);
                end

                Torth=nan(size(XCurr,1),sum(previousBlocks));
                spanStart=[0;cumsum(previousBlocks)]+1;
                spanEnd=cumsum(previousBlocks);
                for j=1:i-1
                    nanMask=~any(ismissing(mdl(j).predict.T),2);
                    Torth(nanMask,spanStart(j):spanEnd(j))=mdl(j).predict.T(nanMask,:);
                end
                nanMask=~any(ismissing(Torth),2);
                XCurr(nanMask,:)=XCurr(nanMask,:)-Torth(nanMask,:)...
                    *pinv(Torth(nanMask,:)'*Torth(nanMask,:))*Torth(nanMask,:)'*XCurr(nanMask,:);
            end
            mdl(i)=DiPLSpredict(XCurr,mdl(i),'Y',Y);%mdl(i)=DiPLSpredict(XCurr(nanMask,:),mdl(i),'Y',Y(nanMask,:));
            predictY(:,:,i)=mdl(i).predict.Y;
            U_predict(size(U_predict,1)-size(mdl(i).predict.U,1)+1:end,1:size(mdl(i).predict.U,2),i)=mdl(i).predict.U;
            C(1:size(mdl(i).Q,1),1:size(mdl(i).Q,2),i)=mdl(i).Q;

        end
        Ypredict=sum(predictY,3);
    otherwise

end
if all(Y==0)
    predict.Y=Ypredict;
else
    nanMask=~any(all(ismissing(U_predict),2),3);
    TSS=    mean((Y(nanMask,:)).^2,1,'omitnan');%-mean(Y(nanMask,:))
    yTemp=  U_predict(nanMask,:,:).*C;
    
    PRESS=  mean((repmat(Y(nanMask,:),[1,size(U_predict,2),size(U_predict,3)])-...
        yTemp).^2,1,'omitnan');
    Q2Y=    1-(PRESS/TSS);
    R2Y=    nan(1,size(Y,2));
    RMSE=   nan(1,size(Y,2));
    MAE=    nan(1,size(Y,2));
    predict.Y=Ypredict;
    nanMask=~isnan(Ypredict);
    for j=1:size(Y,2)
        ESS= mean((Y(:,j)-Ypredict(:,j)).^2,1,'omitnan');
        R2Y(j)= 1-ESS/TSS;
        %         R2Y(j)=corr(Y(:,j),Ypredict(:,j),'Rows','complete');
        RMSE(j)=sqrt(mean((Y(nanMask,j)-Ypredict(nanMask,j)).^2));
        MAE(j)= mean(abs(Y(:,j)-Ypredict(:,j)),'omitnan');
    end
    predict.Q2Y=Q2Y;
    predict.R2Y=R2Y;%.^2;
    predict.Ytrue=Y;
    predict.RMSE=RMSE;
    predict.MAE=MAE;
end
mdl(i+1).predict=predict;