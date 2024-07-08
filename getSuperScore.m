function mdl=getSuperScore(mdl,blocks,target)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                       Superscore calculation
%       WESTERHUIS, Johan A.; KOURTI, Theodora; MACGREGOR, John F.
%       Analysis of multiblock and hierarchical PCA and PLS models.
%     Journal of Chemometrics: A Journal of the Chemometrics Society,
%                        1998, 12.5: 301-321.
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
% Explain here
%-----------------------------------------------------------------<\header>
arguments
    mdl           (1,:) struct   {mustBeUnderlyingType(mdl,'struct')}
    blocks        (1,:) cell     {}
    target        (1,1) string   {} = "base"
end
%define blocks
nBlocks=length(blocks);
switch target
    case "base"
        XBlock=mdl.Xiter(:,:,1);
%         YBlock=mdl.Yiter(:,:,1);
    case "predict"
        XBlock=mdl.predict.X;
%         YBlock=mdl.predict.Y;
    otherwise
        keyboard
end
tss=@(x)sum(sum([(x),-repmat(mean(x,'omitnan'),size(x,1),1)],2,'omitnan').^2,'all');
Rx=zeros(size(mdl.T,2),nBlocks);
tssB=zeros(1,nBlocks);
for k=1:nBlocks
tssB(k)=tss(XBlock(:,blocks{k}(1):blocks{k}(2)));
end
XTemp=XBlock;
blockStep=0:nBlocks:nBlocks*size(mdl.T,2);
Pb=zeros(size(mdl.P));
varPortion=zeros(size(mdl.T,2),nBlocks);
Wb=zeros(size(mdl.W));
Tb=zeros(size(mdl.T,1),size(mdl.T,2)*2);
Wt=zeros(nBlocks,size(mdl.T,2));
nWt=zeros(nBlocks,size(mdl.T,2));
oneWt=zeros(nBlocks,size(mdl.T,2));
nanMask=~any(ismissing(mdl.U),2)&~any(ismissing(mdl.T),2);
for i=1:size(mdl.T,2)
    for j=1:nBlocks
        span = blocks{j}(1):blocks{j}(2);
        Pb(span,i) = XTemp(nanMask,span)'*mdl.T(nanMask,i)/(mdl.T(nanMask,i)'*mdl.T(nanMask,i));
        Wb(span,i) = XTemp(nanMask,span)'*mdl.U(nanMask,i)...
            /(mdl.U(nanMask,i)'*mdl.U(nanMask,i));
        Tb(nanMask,blockStep(i)+j) = XTemp(nanMask,span)*Wb(span,i);
        Eb=XTemp(nanMask,span) - mdl.T(nanMask,i)*Pb(span,i)';
        Eb2=XBlock(nanMask,span) - mdl.T(nanMask,i)*Pb(span,i)';
        varPortion(i,j)=1-(trace(Eb2*Eb2')/trace(XBlock(:,span)*XBlock(:,span)'));
        XTemp(nanMask,span) = Eb;
        Rx(i,j) = (tssB(j) - tss(XTemp(nanMask,span)))/tssB(j);
        
    end
    
    Wt(:,i) = Tb(nanMask,blockStep(i)+1:blockStep(i+1))'*mdl.U(nanMask,i)/(mdl.U(nanMask,i)'*mdl.U(nanMask,i));
    nWt(:,i) = normalize(Wt(:,i),'norm');
    oneWt(:,i) = nWt(:,i)./sum(nWt(:,i));
end

mdl.multiblockStats.Tb=Tb;
mdl.multiblockStats.Pb=Pb;
mdl.multiblockStats.Wb=Wb;
mdl.multiblockStats.Wt=Wt;
mdl.multiblockStats.nWt=nWt;
mdl.multiblockStats.oneWt=oneWt;
mdl.multiblockStats.bRy=(mdl.Ry./sum(mdl.Ry)).*nWt.^2;
mdl.multiblockStats.tbRy=sum(mdl.multiblockStats.bRy,2);
end


