function mdl=SOPLS(X,Y,nComp,algorithm,options)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                               SOPLS
%                     SOPLS algorithm reported in:
%           http://dx.doi.org/
%                  Paper authors: Yining Dong and S. Joe Qin
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% SOPLS model building
% - X matrix
% X matrix to analyze
% - Y vector
% Y vector to analyze
% %%TODO
%-----------------------------------------------------------------<\header>
arguments
    X                      (1,:) cell      {mustBeA(X,"cell"),mustBeNonempty}
    Y                      (:,:) double    {mustBeReal,mustBeNonempty}
    nComp                  (1,:) cell      {mustBeEqualSize(X,nComp)}
    algorithm              (1,:) char      {mustBeMember(algorithm,{'nipals','DiPLS'})}
    options.toll           (1,1) double    {mustBePositive,mustBeReal} =10^-5
    options.maxIter        (1,1) double    {mustBePositive,mustBeReal} =100
    options.isAugmented    (1,1) logical   =true
    options.dynamicOrder   (1,:) double    {mustBeEqualSize(options.dynamicOrder,X)} = ones(size(X))
    options.needsPreprocess(1,1) logical   =true
    options.preprocessType (1,:) char      {mustBeMember(options.preprocessType,{'zscore','block','both','center'})}='zscore';
end
%Blocks in X are already ordered, nComp is a cell array with one value of
%"number of components" per block
%init
mdl(size(X,2))=mdlInit;
switch algorithm
    case 'nipals'
        Yiter=Y;
        for i=1:size(X,2)
            currentBlock=i;
            XCurr=X{currentBlock};
            if i>1
                Torth=nan(size(XCurr,1),sum(cat(1,nComp{1:i-1})));
                spanStart=[0;cumsum(cat(1,nComp{:,:}))]+1;
                spanEnd=cumsum(cat(1,nComp{:,:}));
                
                for j=1:i-1
                    
                    Torth(:,spanStart(j):spanEnd(j))=mdl(j).T;
                end
                
                XCurr=XCurr-Torth*pinv(Torth'*Torth)*Torth'*XCurr;
                % ASDASDASD
            end
            mdl(i)=custom_pls(XCurr,Yiter,nComp{currentBlock},...
                "isFirstblock",i==1,"needsPreprocess",options.needsPreprocess,"preprocessType",options.preprocessType);
            Yiter=mdl(i).Y_i(:,:,end);
        end
        
    case 'DiPLS'
        Yiter=Y;
        for i=1:size(X,2)
            currentBlock=i;
            XCurr=X{currentBlock};
            if i>1
                Xprev=X{currentBlock-1};
                Torth=nan(size(XCurr,1),min(sum(cat(1,nComp{1:i-1})),size(Xprev,2)));
                spanStart=[0;cumsum(min([cat(1,nComp{:,:}),...
                    [size(Xprev,2);size(Xprev,2)]],[],2))]+1;
                spanEnd=cumsum(min([cat(1,nComp{:,:}),...
                    [size(Xprev,2);size(Xprev,2)]],[],2));
                
                for j=1:i-1
                    
                    Torth(:,spanStart(j):spanEnd(j))=mdl(j).T;
                end
                XCurr=XCurr-Torth*pinv(Torth'*Torth)*Torth'*XCurr;
            end
            
            mdl(i)=DiPLS(XCurr,Yiter,nComp{currentBlock},options.dynamicOrder(currentBlock),...
                "toll",options.toll,"maxIter",options.maxIter,"isAugmented",options.isAugmented,...
                "isFirstblock",i==1,"needsPreprocess",options.needsPreprocess,"preprocessType",options.preprocessType);
            Yiter=mdl(i).Y_i(:,:,end);
        end
    otherwise
end





end










% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end
end