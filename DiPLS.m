function mdl=DiPLS(X,Y,nComp,s,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $SODA
%                               DiPLS
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
% this function calculates the parameters for one variable DiPLS X
% - X matrix
% X matrix to analyze
% - Y matrix
% Y matrix to analyze
% - n
% Number of components to extract
% - s
% Autoregressive order of the model
% - "toll",num NAME PAIR ARGUMENT maximum number of iterations
% NAME-VALUE PAIR convergence tollerance
% - "maxIter",num NAME PAIR ARGUMENT maximum number of iterations
% NAME-VALUE PAIR maximum number of iterations
% - "isAugmented",bool
% NAME-VALUE PAIR is the matrix augmented?
%-----------------------------------------------------------------<\header>
defaultToll = 10^-8;
defaultMaxiterations=1000;
defaultisAugmented=false;
defaultneedPreprocess=true;
defaultpreprocessType='zscore';
defaultisFirstblock=true;
checkInput=inputParser;
addRequired(checkInput, 'X')
addRequired(checkInput, 'Y',@(x)mustBeEqualSize(x,X,1))
addRequired(checkInput, 'nComp')
addRequired(checkInput, 's')
addParameter(checkInput,'toll',defaultToll,@isnumeric);
addParameter(checkInput,'maxIter',defaultMaxiterations,@isnumeric)
addParameter(checkInput,'isAugmented',defaultisAugmented,@islogical)
addParameter(checkInput,'needsPreprocess',defaultneedPreprocess,@islogical)
addParameter(checkInput,'preprocessType',defaultpreprocessType,@(x)mustBeMember(x,{'zscore','block','both','center'}))
addParameter(checkInput,'isfirstblock',defaultisFirstblock,@islogical)
parse(checkInput,X,Y,nComp,s,varargin{:})
%init
X=checkInput.Results.X;
Y=checkInput.Results.Y;
nComp=checkInput.Results.nComp;
A=checkInput.Results.s(1);
Ay=checkInput.Results.s(end);
toll=checkInput.Results.toll;
maxIter=checkInput.Results.maxIter;
isAugmented=checkInput.Results.isAugmented;
needsPreprocess=checkInput.Results.needsPreprocess;
preprocessType=checkInput.Results.preprocessType;
isFirstblock=checkInput.Results.isfirstblock;

%--------------------------------------------------------------------------
[~,m] = size (X);           % n = sample size; m = number of input variables
nComp=min(nComp,m);               % maximum number of latent variables
[n,k] = size(Y);            % n=sample size; k = number of output variables
[~,m]=size(X);
if isAugmented
    m=m/(A+1);
    Xorig=X(:,1:m);
else
    Xorig=X;
end
%--------------------------------------------------------Data pre-treatment
% Autoscaling (X and Y)
if needsPreprocess
    switch preprocessType
        case 'center'
            if istable(X)||istimetable(X)
                [X,xCenter]=normalize(X.Variables,'center');
                xScale=1;
            else
                [X,xCenter]=normalize(X,'center');
                xScale=1;
            end
            if isFirstblock
                if istable(Y)||istimetable(Y)
                    [Y,yCenter]=normalize(Y.Variables,'center');
                    yScale=1;
                else
                    [Y,yCenter]=normalize(Y,'center');
                    yScale=1;
                end
            else
                if istable(Y)||istimetable(Y)
                    [~,yCenter]=normalize(Y.Variables,'center');
                    yScale=1;
                else
                    [~,yCenter]=normalize(Y,'center');
                    yScale=1;
                end
            end
        case 'zscore'
            if istable(X)||istimetable(X)
                [X,xCenter,xScale]=normalize(X.Variables,'zscore');
            else
                [X,xCenter,xScale]=normalize(X,'zscore');
            end
            if isFirstblock
                if istable(Y)||istimetable(Y)
                    [Y,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [Y,yCenter,yScale]=normalize(Y,'zscore');
                end
            else
                if istable(Y)||istimetable(Y)
                    [~,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [~,yCenter,yScale]=normalize(Y,'zscore');
                end
            end
        case 'block'
            if istable(X)||istimetable(X)
                xCenter=mean(X.Variables);
                xScale=sum(sum(X.Variables.^2));
                X=(X.Variables-xCenter)./xScale;
            else
                xCenter=mean(X);
                xScale=sum(sum(X.^2));
                X=(X-xCenter)./xScale;
            end
            if isFirstblock
                if istable(Y)||istimetable(Y)
                    [Y,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [Y,yCenter,yScale]=normalize(Y,'zscore');
                end
            else
                if istable(Y)||istimetable(Y)
                    [~,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [~,yCenter,yScale]=normalize(Y,'zscore');
                end
            end
        case 'both'
            if istable(X)||istimetable(X)
                [X,xCenter,xScale]=normalize(X.Variables,'zscore');
                X=X./sum(sum(X.^2));
            else
                [X,xCenter,xScale]=normalize(X,'zscore');
                X=X./sum(sum(X.^2));
            end
            if isFirstblock
                if istable(Y)||istimetable(Y)
                    [Y,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [Y,yCenter,yScale]=normalize(Y,'zscore');
                end
            else
                if istable(Y)||istimetable(Y)
                    [~,yCenter,yScale]=normalize(Y.Variables,'zscore');
                else
                    [~,yCenter,yScale]=normalize(Y,'zscore');
                end
            end
        otherwise
            disp('????')
            keyboard
    end
else
    xCenter=zeros(1,m);
    xScale=ones(1,m);
    yCenter=zeros(1,1);
    yScale=ones(1,1);
end
P=nan(m,nComp);
b=nan(m,k);
Q=nan(k,nComp);
U=nan(n,nComp);
T=nan(n,nComp);
T_s=nan(n,nComp,A+1);
W=nan(m,nComp);
BetaX=nan(A+1,nComp);
BetaY=nan(Ay+1,nComp);
BetaXhat=nan(A+1,nComp);
Xiter=nan(size(X,1),size(X,2),A+1);
Yiter=nan(size(Y,1),size(Y,2),A+1);
X_i=X;
Y_i=Y;

for i=1:nComp
    Xiter(:,:,i)=X_i;
    Yiter(:,:,i)=Y_i;
    [P(:,i),Q(:,i),U(:,i), T(:,i),W(:,i),BetaX(:,i),BetaY(:,i),BetaXhat(:,i),X_i,Y_i,T_s(:,i,:)] = ...
        dynamicinnernipals(X_i,Y_i,'A',s,'isAugmented',isAugmented,'toll',toll,'maxIter',maxIter);
    
end
%------------------------------Amount of variability explained by the model
Q_calc=     permute(Q,        [3,1,2]);%Q transposed
T_calc=     permute(T,        [1,3,2]);
P_calc=     permute(P,        [3,1,2]);%P transposed
Ts_calc=    permute(T_s,      [1,3,2]);
Beta_calc=  permute(BetaXhat, [1,3,2]);
Xpred=pagemtimes(T_calc,P_calc);
Ypred=pagemtimes(pagemtimes(Ts_calc,Beta_calc),Q_calc);
if isAugmented
    cubeESSx=sum(sum((Xorig-Xpred).^2,[],'omitnan'),2);
    ESSx=permute(cubeESSx,[1,3,2]);
    TSSx=sum(sum((Xorig-mean(Xorig)).^2,1,'omitnan'));
else
    cubeESSx=sum(sum((X-Xpred).^2,[],'omitnan'),2);
    ESSx=permute(cubeESSx,[1,3,2]);
    TSSx=sum(sum((X-mean(X)).^2,1,'omitnan'));
end
ESSy=permute(sum((Y(s+1:end,:)-Ypred(s+1:end,:)).^2,1,'omitnan'),[3,2,1]);
TSSy=sum((Y(s+1:end,:)-mean(Y(s+1:end,:))).^2,1,'omitnan');
Ry=1-(ESSy./TSSy);
Rx=1-(ESSx./TSSx);
RRy=cumsum(Ry);
RRx=cumsum(Rx);
%-------------------------------------------------Build the control indices
mm=size(T,1);
aa=size(T,2);
kk=size(P,1);
Eig=diag(T'*T)/(size(T,1)-1);%diag(diag(((T)'*(T))/(mm-1))+eps);
f=sqrt((1./Eig)+eps);
Ht2=sum((T*diag(f)).^2,2);
Ht2Cont=X*P/diag(power(Eig,0.5))*P';
HotellingLim95F=  (aa*(mm-1)/(mm-aa))*finv(0.95,aa,mm-aa);
HotellingLim95B=((mm-1)^2/mm)*betainv(0.95,(aa/2),.5*(mm-aa-1));
SPE=diag(X*(eye(kk)-P*P')*X').^2;%diag(X*(eye(kk)-P*P')*X');
SPELim95=QLim(SPE,0.95);
if mm<50
    HotellingLim95=HotellingLim95B;
else
    HotellingLim95=HotellingLim95F;
end
%-----------------------------------------------Build the reponse structure
for i=1:size(Y,2)
    bTemp=P'*W;
    bTemp(bTemp==0)=eps;
    bCoeff= cumsum(W*(bTemp.\diag(Q(i,:)')),2);
    b(:,i)= bCoeff(:,end);
end
mdl=mdlInit;
mdl.P               =P;
mdl.Q               =Q;
mdl.U               =U;
mdl.T               =T;
mdl.W               =W;
mdl.C               =[];
mdl.B               =[];
mdl.BetaX           =BetaX;
mdl.BetaY           =BetaY;
mdl.BetaXhat        =BetaXhat;
mdl.X_i             =X_i;
mdl.Xiter           =Xiter;
mdl.Yiter           =Yiter;
mdl.Y_i             =Y_i;
mdl.T_s             =T_s;
mdl.Ry              =Ry;
mdl.Rx              =Rx;
mdl.RRy             =RRy;
mdl.RRx             =RRx;
mdl.bCoeff          =b;
mdl.xScale          =xScale;
mdl.xCenter         =xCenter;
mdl.yScale          =yScale;
mdl.yCenter         =yCenter;
mdl.isAugmented     =isAugmented;
mdl.Hotelling       =Ht2;
mdl.HotellingC      =Ht2Cont;
mdl.HotellingLim95  =HotellingLim95;
mdl.SPE             =SPE;
mdl.SPELim95        =SPELim95;
mdl.predict         =[];
mdl.type            ="DiPLS";
end

%--------------------------------------------------------------------------
%   Helper functions
%--------------------------------------------------------------------------
% Custom validation function
function mustBeEqualSize(a,b,d)
% Test for equal size
if ~isequal(size(a,d),size(b,d))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end
end

%--------------------------------------------------------------------------
function SPElim=QLim(SPE,cl)
theta1 = trace(SPE*SPE'/(size(SPE,1)-1));
theta2 = trace((SPE*SPE'/(size(SPE,1)-1)).^2);
theta3 = trace((SPE*SPE'/(size(SPE,1)-1)).^3);
h0     = 1-2*theta1*theta3/3/(theta2.^2);
if h0<0.001
    h0 = 0.001;
end
ca    = sqrt(2)*erfinv(2*(1-cl));
h1    = ca*sqrt(2*theta2*h0.^2)/theta1;
h2    = theta2*h0*(h0-1)/(theta1.^2);
SPElim = theta1*(1+h1+h2).^(1/h0);
end
