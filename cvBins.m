function [blockSet,permutationSet]=cvBins(x,n,varargin)
defaultReduced=true;
checkInput=inputParser;
addRequired(checkInput,'x',@(x)isnumeric(x)&&ismatrix(x))
addRequired(checkInput,'n',@(x)isnumeric(x)&&isscalar(x))
addParameter(checkInput,'Reduced',defaultReduced,@islogical)
parse(checkInput,x,n,varargin{:})
%init
x=checkInput.Results.x;
n=checkInput.Results.n;
Reduced=checkInput.Results.Reduced;
bin=discretize(1:size(x,1),linspace(0,size(x,1)+1,n+1),'IncludedEdge','right');
blockSet=cell(1,n);
for i=1:n
blockSet{i}=find(bin==i)';
end
if Reduced
permutationSet=nchoosek(2:n-1,ceil((n-2)/2));
else
permutationSet=nchoosek(1:n,ceil(n/2));
end