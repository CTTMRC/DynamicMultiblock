function mdl=SODIPLS(X,Y,nComp,S,options)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $SODA
%                               SO-DiPLS
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
% X matrix to analyze, must be a cell vector with different blocks in
% different cells. Must not be augumented, put only T0 here.
% - Y matrix
% Y matrix to analyze
% -nComp
%Number of components, needs to be a cell; same size as X. 
% - S
% Dynamic order of the model
%-----------------------------------------------------------------<\header>
arguments
    X                      (1,:) cell      {mustBeA(X,"cell"),mustBeNonempty}
    Y                      (:,:) double    {mustBeReal,mustBeNonempty}
    nComp                  (1,:) cell      {mustBeA(nComp,"cell"),mustBeNonempty}
    S                      (1,:) double    {mustBeReal,mustBeNonempty}
    options.toll           (1,1) double    {mustBePositive,mustBeReal} =10^-5
    options.maxIter        (1,1) double    {mustBePositive,mustBeReal} =100
    options.needsPreprocess(1,1) logical   =false
    options.preprocessType (1,:) char      {mustBeMember(options.preprocessType,{'zscore','block','both'})}='zscore';
    options.calibrate      (1,1) logical   =true;
end
%<body>
mdl= SOPLS(X,Y,nComp,'DiPLS','dynamicOrder',S,'isAugmented',false,'needsPreprocess',options.needsPreprocess,'preprocessType',options.preprocessType);

%</body>
end
