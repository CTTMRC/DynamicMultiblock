function [p,q,u,t,w,arX,arY,b,E_x,E_y,t_s] = dynamicinnernipals (X,Y,varargin)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                           Dyamic Inner NIPALS
%                     DiNIPALS algorithm reported in:
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
%  -"A",num
% NAME-VALUE PAIR Autoregressive order of the model (valid for x)
% - "toll",num NAME PAIR ARGUMENT maximum number of iterations
% NAME-VALUE PAIR convergence tollerance
% - "maxIter",num NAME PAIR ARGUMENT maximum number of iterations
% NAME-VALUE PAIR maximum number of iterations
% - "isAugmented",bool
% NAME-VALUE PAIR is the matrix augmented?
%-----------------------------------------------------------------<\header>
defaultToll = 10^-5;
defaultAutoregressiveorder = 1;
defaultMaxiterations=100;
defaultStartingARX=[];
defaultisAugmented=true;
checkInput=inputParser;
addRequired(checkInput,'X')
addRequired(checkInput,'Y')
addParameter(checkInput,'A',defaultAutoregressiveorder,@isnumeric);
addParameter(checkInput,'toll',defaultToll,@isnumeric);
addParameter(checkInput,'maxIter',defaultMaxiterations,@isnumeric)
addParameter(checkInput,'isAugmented',defaultisAugmented,@islogical)
addParameter(checkInput,'startingBeta',defaultStartingARX,@isnumeric)
parse(checkInput,X,Y,varargin{:})
X=checkInput.Results.X;
Y=checkInput.Results.Y;
A=checkInput.Results.A;
toll=checkInput.Results.toll;
maxIter=checkInput.Results.maxIter;
isAugmented=checkInput.Results.isAugmented;
arX=checkInput.Results.startingBeta;
%Initialization
[~,m]=size(X);
if isAugmented
    m=m/(A+1);
    Xorig=X(:,1:m);
else
    Xorig=X;
    [X,Y]=matrixAugment(X,Y,A);
end
Y(isnan(Y),:)=0;
i_0=randi(width(Y));
u = Y(:,i_0);
w=ones(m,1);
w_old=zeros(m,1);
if isempty(arX)
    arX=[1;zeros(A,1)];
end
arY=[1;zeros(A,1)];
currIter=0;
% V=zeros(maxIter,1); %to follow convergence
% Conv=zeros(maxIter,1); %to follow convergence
while norm(abs(w)-abs(w_old))>toll && currIter<=maxIter    %convergence check
    currIter=currIter+1;
    w_old=w;
    w=(kron(arX,eye(m)))'*X'*u;%w=(Beta⊗I)'Z'Yq=(Beta⊗I)'Z'u
    w=w/norm(w);
    t = X*kron(eye(A+1),w);%t=Z(I⊗w)
    q=Y'*X*kron(arX,w);%q=Y'Z(Beta⊗w);
    q=q/norm(q);
    u=Y*q;
    arX=t'*u;%Beta=(I⊗w)'Z'Yq =(I⊗w)'Z'u = t'u;
    arX=arX/norm(arX);
    %%%%%HERE FOR CHECKING ITERATION AND CONVERGENCE%%%%%%%%%%%%%%%%%%%%%%%
    %     V(currIter,1)=q'*Y'*X*kron(arX,w); %to follow convergence
    %     Conv(currIter,1)=norm(w-w_old); %to follow convergence
    %     fprintf(' targetFun outer:%f;\n targetFun inner:%f;\n',...
    %         V(currIter,1),Conv(currIter,1)); %to follow convergence
    %     fprintf(' arO:%f\n',arX); %to follow convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
q = Y'*u/(u'*u); q=q/norm(q);           % Y-loading
u = t*((t'*t)\t'*Y*q);                  % U_hat
b = (t'*t)\t'*u; %b=b/norm(b);          % BetaHat
t_s=t;                                  % Autoregressive X-score
t = Xorig*w;                            % X-score
p = Xorig'*t/(t'*t);                    % X-loading
E_x = Xorig-t*p';                       % Residual matrix for X
E_y = Y-u*q';                           % Residual matrix for Y

if isAugmented
    E_x=X-t_s*kron(p,eye(A+1))';
else
    E_y=[nan(A,size(Y,2));E_y];
    t_s=[nan(A,A+1);t_s];
    u=[nan(A,1);u];
end
end
%--------------------------------------------------------------------------
%   Helper functions
%--------------------------------------------------------------------------
function [augmentedX,augmentedY]=matrixAugment(X,Y,ArO)
%Augmented Matrix needs to be X(t),X(t-1),...,X(t-ArO) for algorithmic
%purposes
[~,m]=size(X);
augmentedX=zeros(size(X,1)-ArO,size(X,2)*(ArO+1));
augmentedY=Y(ArO+1:end,:);
augmentedX(:,1:size(X,2))=X(ArO+1:end,:);
for i=1:ArO
    startIdx=(m*i)+1;
    endIdx=m*(i+1);
    augmentedX(:,startIdx:endIdx)=X(ArO+1-i:end-i,:);
end
end