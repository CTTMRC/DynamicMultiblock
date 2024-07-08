function mdl=SODA(X,Y,nComp,options)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $SODA
%                               SODA-PLS
%                     SODA-PLS algorithm reported in:
%                               DOI-HERE
%                  Paper authors: M.Cattaldo, A.Ferrer, I.Måge
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo, Ingrid Måge
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% this function calculates the parameters for one variable DiPLS X
% - X matrix
% X matrix to analyze, must be a cell vector with different blocks in
% different cells
% - Y matrix
% Y matrix to analyze
% -nComp
%Number of components, needs to be a cell; same size as X. 
%-----------------------------------------------------------------<\header>
arguments
    X                      (1,:) cell      {mustBeA(X,"cell"),mustBeNonempty}
    Y                      (:,:) double    {mustBeReal,mustBeNonempty}
    nComp                   (1,:) cell      {mustBeA(nComp,"cell"),mustBeNonempty}
    options.toll           (1,1) double    {mustBePositive,mustBeReal} =10^-5
    options.maxIter        (1,1) double    {mustBePositive,mustBeReal} =100
    options.needsPreprocess(1,1) logical   =false
    options.preprocessType (1,:) char      {mustBeMember(options.preprocessType,{'zscore','block','both'})}='zscore';
    options.calibrate      (1,1) logical   =true;
end
%<body>
mdl= SOPLS(X,Y,nComp,'nipals',"needsPreprocess",options.needsPreprocess,"preprocessType",options.preprocessType);


%</body>
end
% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eID = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eID,msg))
end
end