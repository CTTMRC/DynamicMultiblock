function mdl = custom_pls(X,Y,nComp,options)
%------------------------------------------------------------maca/gabrieleb
%------------------------------------------------------------------<header>
%                      Partial Least Squares (PLS)
%                         NIPALS algorithm
%                         In-house function
%                           Gabriele Bano
%         CAPE-Lab - Computer Aided Process Engineering Laboratory
%                      University of Padova (Italy)
%                     Modified by: Marco Cattaldo
%                           Nofima (norway)
%-----------------------------------------------------------------<\header>
arguments
    X                      (:,:) double    {mustBeReal,mustBeNonempty}
    Y                      (:,:) double    {mustBeReal,mustBeNonempty}
    nComp                  (1,1) double    {mustBeReal,mustBeNonempty}
    options.toll           (1,1) double    {mustBePositive,mustBeReal} =10^-5
    options.maxIter        (1,1) double    {mustBePositive,mustBeReal} =100
    options.needsPreprocess(1,1) logical   =false
    options.preprocessType (1,:) char      {mustBeMember(options.preprocessType,{'zscore','block','both','center'})}='zscore';
    options.isFirstblock   (1,1) logical   =true
end
%--------------------------------------------------------------------------
[~,m] = size (X);         % n = sample size; m = number of input variables
[n,k] = size(Y);          % n=sample size; k = number of output variables
%--------------------------------------------------------Data pre-treatment
% Mean centering and autoscaling (X and Y)
% Autoscaling (X and Y)
needsPreprocess=options.needsPreprocess;
preprocessType=options.preprocessType;
isFirstblock=options.isFirstblock;
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
%--------------------------------------------------------------------------
% X- and Y- scores/loadings computation (NIPALS algorithm)
%--------------------------------------------------------------------------
%inizialitation
toll    = options.toll ;  %tolerance on residuals
maxIter =options.maxIter;
P       =nan(m,nComp);
Q       =nan(k,nComp);
U       =nan(n,nComp);
T       =nan(n,nComp);
W       =nan(m,nComp);
C       =nan(k,nComp);
B       =nan(1,nComp);
b       =nan(m,k);
X_i     =nan(n,m,nComp);
Y_i     =nan(n,k,nComp);
Xiter   =nan(n,m,nComp);
Yiter   =nan(n,k,nComp);
for j = 1 : nComp
    if j==1
        Xiter(:,:,j) = X;
        Yiter(:,:,j) = Y;
    else
        Xiter(:,:,j) = X_i(:,:,j-1);
        Yiter(:,:,j) = Y_i(:,:,j-1);
    end
    [P(:,j),Q(:,j),U(:,j), T(:,j),W(:,j),C(:,j),B(:,j),X_i(:,:,j),Y_i(:,:,j)] = nipals (Xiter(:,:,j),Yiter(:,:,j),toll,maxIter);
    
end

%------------------------------Amount of variability explained by the model
C_calc=permute(C,[3,1,2]);
T_calc=permute(T,[1,3,2]);
P_calc=permute(P,[3,1,2]);
B_calc=permute(B,[1,3,2]);
Xpred=pagemtimes(T_calc,P_calc);
Ypred=pagemtimes(B_calc,pagemtimes(T_calc,C_calc));
ESSy=permute(sum(sum((Y-Ypred).^2,[],'omitnan'),1),[1,3,2]);
ESSx=permute(sum(sum((X-Xpred).^2,[],'omitnan'),2),[1,3,2]);
TSSy=permute(sum((Y-mean(Y)).^2,1,'omitnan'),[1,3,2]);
TSSx=sum(sum((X-mean(X)).^2,1,'omitnan'));
Ry=1-(pagemtimes(ESSy,(1./TSSy)));
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
%-------------------------------------------------------------------------
% Final dimensions of the matrices computed in the loop
%
% W: dimension m*max_lvs
% T: dimension n*max_lvs
% C: dimension k*max_lvs
% U: dimension n*max_lvs
% P: dimension m*max_lvs
% Q: dimension k*max_lvs
% B: dimension 1*max_lvs
%--------------------------------------------------------------------------
for i=1:size(Y,2)
    bTemp=P'*W;
    bTemp(bTemp==0)=eps;
    bCoeff= W/bTemp*B';
    b(:,i)= bCoeff(:,end);
end
mdl=mdlInit;
mdl.P               =P;
mdl.Q               =Q;
mdl.U               =U;
mdl.T               =T;
mdl.W               =W;
mdl.C               =C;
mdl.B               =B;
mdl.X_i             =X_i;
mdl.Y_i             =Y_i;
mdl.Xiter           =Xiter;
mdl.Yiter           =Yiter;
mdl.xScale          =xScale;
mdl.xCenter         =xCenter;
mdl.yScale          =yScale;
mdl.yCenter         =yCenter;
mdl.Ry              =Ry;
mdl.Rx              =Rx;
mdl.RRy             =RRy;
mdl.RRx             =RRx;
mdl.bCoeff          =b;
mdl.BetaX           =[];
mdl.BetaY           =[];
mdl.BetaXhat        =[];
mdl.T_s             =[];
mdl.isAugmented     =[];
mdl.Hotelling       =Ht2;
mdl.HotellingC      =Ht2Cont;
mdl.HotellingLim95  =HotellingLim95;
mdl.SPE             =SPE;
mdl.SPELim95        =SPELim95;
mdl.predict         =[];
mdl.type            ="nipals";
end
%-----------------------------------------------------------------------------------------------
% NIPAlS Algorithm (Ref. A. Hoskuldsson- PLS regression method - Journal of Chemometrics (1988))
%-----------------------------------------------------------------------------------------------
function [p,q,u,t,w,c,b,E_x,E_y] = nipals (X,Y,toll,n)
%Initialization
nanMask=~any(isnan(X),2)&~any(isnan(Y),2);
XOrig=X;
YOrig=Y;
X=X(nanMask,:);
Y=Y(nanMask,:);
u = ones(size(Y(:,1)));       %initial value for u just to enter in the first cycle)
u_new = Y(:,1);
iter=0;
while norm(u_new-u)>toll && iter<n  %convergence check
    u = u_new;
    w = X'*u/(u'*u);                % weight    : dimension : m*1
    w = w/norm(w);
    t = X*w;                        % x-loading : dimension : n*1
    c = Y'*t/(t'*t);                %           : dimension : k*1
    c = c/norm(c);
    u_new = Y*c/(c'*c);             % y-score   : dimension : n*1
    iter=iter+1;
end
u   = u_new;
q   = Y'*u/(u'*u);
p   = X'*t/(t'*t);
b   = u'*t/(t'*t);
E_x = X-t*p';
E_y = Y-b*t*c';
u   = YOrig*c/(c'*c);
t   = XOrig*w;
XOrig(nanMask,:)=E_x;
YOrig(nanMask,:)=E_y;
E_x =XOrig;
E_y =YOrig;
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
