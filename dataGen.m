function [trainSet,yTrain,calSet,yCalibration,testSet,yTest]=dataGen(weightCase,blockWeight,whichSpectra,noiseLevel,options)
%----------------------------------------------------------------------maca
%------------------------------------------------------------------<header>
%Part of $packageName
%                       Syntetic Data Generation
%                         In-house function
%                           version: v.01
%                       ogAuthor: Marco Cattaldo
%                   lastModified: Marco Cattaldo
%                               NOFIMA
%--------------------------------------------------------------------------
%-------------------------------------------------------------------explain
%The first number refers to the Process block, the second number refers to
%the spectra block
%blockWeight    - must be between 1 and 3: weight of the Contribution
%   1=[0.5 0.5], 2=[0.2 0.8], 3=[0.8 0.4]
%weightCase     - must be between 1 and 4: strength of the Dynamics
%   1=[0,0], 2=[1,1], 3=[2,0], 4=[0,2]
%whichSpectra   - must be between 1 and 2: which spectra scores are used to
%calculate the spectra contribution
%   1=PC1, 2=PC[2-5]
%noiseLevel     - must be between 1 and 2: how much nois there is - more or
%less
%   1=5%, 2=15%
%-----------------------------------------------------------------<\header>
arguments
    weightCase           (1,1) double      {mustBeGreaterThanOrEqual(weightCase,1),mustBeLessThanOrEqual(weightCase,4)} = 1
    blockWeight          (1,1) double      {mustBeGreaterThanOrEqual(blockWeight,1),mustBeLessThanOrEqual(blockWeight,3)} = 1
    whichSpectra         (1,1) double      {mustBeGreaterThanOrEqual(whichSpectra,1),mustBeLessThanOrEqual(whichSpectra,2)} = 1
    noiseLevel           (1,1) double      {mustBeGreaterThanOrEqual(noiseLevel,1),mustBeLessThanOrEqual(noiseLevel,2)} = 1
    options.blockScaling (1,:)             {mustBeMember(options.blockScaling,...
        {'None','SNN','SBS','HBS','SHBS','SBVS','HBVS','SHBVS','SBRS','HBRS','SHBRS'})}= 'None'
    
end

%% Settings
%%SEED
randomSeed=1;
randomGenerator='twister';
blockScaling=char(options.blockScaling);
%%RELATIVEWEIGHT
% weightCase=dFF(z,1);
% blockWeight=dFF(z,2);
% %%MISCELLANEA
% whichSpectra=dFF(z,3);
% noiseLevel
%% Loading Sets
load('.\Data\NIRloads.mat','NIRloads');
load('.\Data\NIRmean.mat', 'NIRmean' );
%% Build Datasets

switch weightCase
    case 1
        processWeightT0=1       ;
        processWeightT1=0       ;
        processWeightT2=0       ;
        spectraWeightT0=1       ;
        spectraWeightT1=0       ;
        spectraWeightT2=0       ;
    case 2
        processWeightT0=1       ;
        processWeightT1=1       ;
        processWeightT2=0       ;
        spectraWeightT0=1       ;
        spectraWeightT1=1       ;
        spectraWeightT2=0       ;
    case 3
        processWeightT0=1       ;
        processWeightT1=1       ;
        processWeightT2=1       ;
        spectraWeightT0=1       ;
        spectraWeightT1=0       ;
        spectraWeightT2=0       ;
    case 4
        processWeightT0=1       ;
        processWeightT1=0       ;
        processWeightT2=0       ;
        spectraWeightT0=1       ;
        spectraWeightT1=1       ;
        spectraWeightT2=1       ;
    otherwise
        keyboard
        
end
%{
% nirMask=~any(ismissing(BC22XNIR(:,1:141)),2);
% trueSpectraBlock=BC22XNIR(nirMask,1:141);
% trueSpectraBlock{:,:}=snv(trueSpectraBlock{:,:});
% trueSpectraBlock=normalize(trueSpectraBlock,'center');
% pcaOpt.display='none';
% pcaOpt.plots='none';
% pcaOpt.preprocessing='';
% pcaNIR=pca(trueSpectraBlock.Variables,5,pcaOpt);
% NIRloads=pcaNIR.loads{2};
%}
P =[0.8691,0.4306,0.2433 ;
    0.3742,0.0468,0.9262 ;
    0.4729,0.3730,0.7983 ;
    0.4010,0.3873,0.8302 ;
    0.4882,0.8287,0.2736];
C11 =[0.5437;
    0.3596;
    0.5342;
    0.3458;
    0.4125];
C12=[0.6858;
    0.2658;
    0.3490;
    0.3977;
    0.4233];
C13=[0.2828;
    0.0698;
    0.6490;
    0.1977;
    0.2233];
A11= [0.6767,   0.5809,    0.9315;
    1.2812,  -0.5343,   -1.6000;
    -1.5083,   0.9991,    0.7529];
A12= [0.7155,  -0.0652,    1.1192;
    1.1132,  -0.5371,   -0.1691;
    - 0.5571,  -1.0748,    0.2330];
A21=[ 0.8528,  -0.3192,   -0.7924,   -0.9598,    0.2987;
    -0.2016,   0.7423,    0.1786,   -0.2509,   -0.8902;
    0.7378,  -1.1877,    0.3143,   -0.5012,    1.8799;
    -1.6747,   1.8337,   -0.1980,    0.5260,   -0.6289;
    0.5729,  -0.7997,    0.1624,    0.0066,    2.5151];
A22=[-0.4931,  -0.8300,   -1.2786,    0.7055,    2.2875;
    0.9962,   0.1601,   -0.8591,    1.9183,    0.0688;
    -1.0858,  -1.1235,    0.5436,   -0.1136,    0.3865;
    -0.9385,   0.8455,    0.1984,    1.1474,   -0.5123;
    0.4630,  -1.0825,   -1.4736,   -0.0482,   -1.6528];

C21=conv(NIRloads(:,1),NIRloads(:,3)./2,'same');
C22=conv(NIRloads(:,2),NIRloads(:,4)./2,'same');
C23=flipud(conv(NIRloads(:,3),NIRloads(:,5)./2,'same'));
C11=(C11./norm(C11)).*processWeightT0;
C12=(C12./norm(C12)).*processWeightT1;
C13=(C13./norm(C13)).*processWeightT2;
C21=(C21./norm(C21)).*spectraWeightT0;
C22=(C22./norm(C22)).*spectraWeightT1;
C23=(C23./norm(C23)).*spectraWeightT2;
%% completely syntetic dataset
rng(randomSeed,randomGenerator);
peStream=normrnd(0     ,0.050  ,1002,5);
seStream=normrnd(0     ,0.0250   ,1002,141);
yStream=normrnd(0      ,1      ,1002,1);
S=[0.1972, 0.0816, 0.0560, 0.0263, 0.0252];
V=NIRloads;
%%%%DYNAMICS ON X
pfStream=normrnd(0     ,0.050  ,1002,3);
sfStream=normrnd(0     ,0.025  ,1002,5);
pStream=normrnd(0      ,1      ,1004,3);
sStream=normrnd(0      ,1      ,1004,5);
pStreamDyn= pStream(2:1003,:)*A11+pStream(3:1004,:)*A12+pfStream;
sStreamDyn= sStream(2:1003,:)*A21+sStream(3:1004,:)*A22+sfStream;
spectraBlock=sStreamDyn.*S*V'+repmat(NIRmean,[1002,1])+seStream;
productBlock=pStreamDyn*P'+peStream;
switch whichSpectra
    case 1
        spectraYBlock=sStreamDyn(:,1).*S(1)*V(:,1)'+repmat(NIRmean,[1002,1]);
    case 2
        spectraYBlock=sStreamDyn(:,2:end).*S(2:end)*V(:,2:end)'+repmat(NIRmean,[1002,1]);
    otherwise
        keyboard
end
%NO DYNAMICS ON X
% pStream=normrnd(0      ,4      ,1002,3);
% sStream=normrnd(0      ,1      ,1002,5);
% spectraBlock=sStream.*S*V'+repmat(mean(BC22XNIR{:,1:141},1,'omitnan'),[1002,1]);
% switch whichSpectra
%     case 1
%     spectraYBlock=sStream(:,1).*S(1)*V(:,1)'+repmat(mean(BC22XNIR{:,1:141},1,'omitnan'),[1002,1]);
%     case 2
%     spectraYBlock=sStream(:,2:end).*S(2:end)*V(:,2:end)'+repmat(mean(BC22XNIR{:,1:141},1,'omitnan'),[1002,1]);
%     otherwise
%         keyboard
% end
% productBlock=pStream*P'+eStream;
t0=struct;
t1=struct;
t2=struct;
t0.processBlock             =productBlock;
t0.spectraBlock             =spectraBlock;
t0.spectraYBlock            =spectraYBlock;
t1.processBlock             =nan(size(productBlock));
t1.spectraBlock             =nan(size(spectraBlock));
t1.spectraYBlock            =nan(size(spectraYBlock));
t2.processBlock             =nan(size(productBlock));
t2.spectraBlock             =nan(size(spectraBlock));
t2.spectraYBlock            =nan(size(spectraYBlock));
t1.processBlock(2:end,:)    =t0.processBlock(1:end-1,:);
t2.processBlock(3:end,:)    =t0.processBlock(1:end-2,:);
t1.spectraBlock(2:end,:)    =t0.spectraBlock(1:end-1,:);
t2.spectraBlock(3:end,:)    =t0.spectraBlock(1:end-2,:);
t1.spectraYBlock(2:end,:)   =t0.spectraYBlock(1:end-1,:);
t2.spectraYBlock(3:end,:)   =t0.spectraYBlock(1:end-2,:);
yProcess=t0.processBlock*C11...
    +t1.processBlock*C12...
    +t2.processBlock*C13;
ySpectra=t0.spectraYBlock*C21...
    +t1.spectraYBlock*C22...
    +t2.spectraYBlock*C23;
tss=@(x)sum(sum([(x),-repmat(mean(x,'omitnan'),size(x,1),1)],2,'omitnan').^2,'all');
fsolveOpt=optimoptions(@fsolve,'Display','none');
weightFun=@(x,a,b,target)((tss(a*(x)))/(tss(b*(1-x))))-target;
optWeight=@(x,target)weightFun(x,yProcess,ySpectra,target);
%super hard block rank scaling
switch blockWeight
    case 1
        case1=@(x)optWeight(x,1);
        weight=fsolve(case1,0.5,fsolveOpt);
        ySynth=  yProcess*(weight)...
            +ySpectra*(1-weight);
    case 2
        case2=@(x)optWeight(x,0.20);
        weight=fsolve(case2,0.5,fsolveOpt);
        ySynth=  yProcess*(weight)...
            +ySpectra*(1-weight);
    case 3
        case3=@(x)optWeight(x,4);
        weight=fsolve(case3,0.5,fsolveOpt);
        ySynth=  yProcess*(weight)...
            +ySpectra*(1-weight);
    otherwise
        keyboard
end

yStd=std(ySynth,[],'omitnan');
EYcorr=corr(yStream(3:end),ySynth(3:end));
fsolveOpt=optimoptions(@fsolve,'Display','none');
switch noiseLevel
    case 1
        f=@(x)1+x+(2*sqrt(x)*EYcorr)-1.0526^2;
        mult=sqrt(fsolve(f,0.05,fsolveOpt));
        ySynth=ySynth+(yStream*mult*yStd);
    case 2
        f=@(x)1+x+(2*sqrt(x)*EYcorr)-1.1579^2;
        mult=sqrt(fsolve(f,0.2,fsolveOpt));
        ySynth=ySynth+(yStream*mult*yStd);
    otherwise
        keyboard
end
sets.Train=3:602;
sets.Calibration=3:602;
sets.Test=603:1002;
%%
[xPtrain0,pCenter,pScale]=  normalize(t0.processBlock(sets.Train,:),'zscore');
[xStrain0,sCenter]=         normalize(t0.spectraBlock(sets.Train,:),'center');
xPtrain1=                   normalize(t1.processBlock(sets.Train,:),'center',pCenter,'scale',pScale);
xStrain1=                   normalize(t1.spectraBlock(sets.Train,:),'center',sCenter,'scale',1);
xPtrain2=                   normalize(t2.processBlock(sets.Train,:),'center',pCenter,'scale',pScale);
xStrain2=                   normalize(t2.spectraBlock(sets.Train,:),'center',sCenter,'scale',1);
[yTrain,yCenter,yScale]=    normalize(ySynth(sets.Train,:),'zscore');
xPcalibration0=             normalize(t0.processBlock(sets.Calibration,:),'center',pCenter,'scale',pScale);
xScalibration0=             normalize(t0.spectraBlock(sets.Calibration,:),'center',sCenter,'scale',1);
xPcalibration1=             normalize(t1.processBlock(sets.Calibration,:),'center',pCenter,'scale',pScale);
xScalibration1=             normalize(t1.spectraBlock(sets.Calibration,:),'center',sCenter,'scale',1);
xPcalibration2=             normalize(t2.processBlock(sets.Calibration,:),'center',pCenter,'scale',pScale);
xScalibration2=             normalize(t2.spectraBlock(sets.Calibration,:),'center',sCenter,'scale',1);
yCalibration=               normalize(ySynth(sets.Calibration,:),'center',yCenter,'scale',yScale);
xPtest0=                    normalize(t0.processBlock(sets.Test,:),'center',pCenter,'scale',pScale);
xStest0=                    normalize(t0.spectraBlock(sets.Test,:),'center',sCenter,'scale',1);
xPtest1=                    normalize(t1.processBlock(sets.Test,:),'center',pCenter,'scale',pScale);
xStest1=                    normalize(t1.spectraBlock(sets.Test,:),'center',sCenter,'scale',1);
xPtest2=                    normalize(t2.processBlock(sets.Test,:),'center',pCenter,'scale',pScale);
xStest2=                    normalize(t2.spectraBlock(sets.Test,:),'center',sCenter,'scale',1);
yTest=                      normalize(ySynth(sets.Test,:),'center',yCenter,'scale',yScale);
% pTSS=                       @(x)sqrt(sum(std(x)));%sum(std(xPtrain0));%sqrt(sum(var(xPtrain0)));%sqrt(size(xPtrain0,2));%sqrt(sum(std(xPtrain0).^2));%norm(xPtrain0);%
% sTSS=                       @(x)sqrt(sum(std(x)));%sum(std(xStrain0));%sqrt(sum(var(xStrain0)));%sqrt(size(xStrain0,2));%sqrt(sum(std(xStrain0).^2));%norm(xStrain0);%
% pTrainset=[{xPtrain0./pTSS(xPtrain0)},{xPtrain1./pTSS(xPtrain0)},{xPtrain2./pTSS(xPtrain0)}];
% sTrainset=[{xStrain0./sTSS(xStrain0)},{xStrain1./sTSS(xStrain0)},{xStrain2./sTSS(xStrain0)}];
% pCalset=[{xPcalibration0./pTSS(xPtrain0)},{xPcalibration1./pTSS(xPtrain0)},{xPcalibration2./pTSS(xPtrain0)}];
% sCalset=[{xScalibration0./sTSS(xStrain0)},{xScalibration1./sTSS(xStrain0)},{xScalibration2./sTSS(xStrain0)}];
% pTestset=[{xPtest0./pTSS(xPtrain0)},{xPtest1./pTSS(xPtrain0)},{xPtest2./pTSS(xPtrain0)}];
% sTestset=[{xStest0./pTSS(xStrain0)},{xStest1./pTSS(xStrain0)},{xStest2./pTSS(xStrain0)}];
switch blockScaling
    % ref listed in Data preprocessing for multiblock modelling –
    % A systematization with new methods DOI:https://doi.org/10.1016/j.chemolab.2020.103959
    case 'None'
        %no block scaling
        N  = @(x)1./ones(1,size(x,2));
        nN = @(x)x./ones(1,size(x,2));
        f  = {N,nN};
    case 'SNN'
        %block scale to unit norm
        SSN   = @(x) 1./sqrt(sum(sum(x.^2)));
        nSSN  = @(x) x./sqrt(sum(sum(x.^2)));
        f     = {SSN,nSSN};
    case 'SBS'
        %soft block scaling
        SBS   = @(x) 1./nthroot(size(x,2),4);
        nSBS  = @(x) x./nthroot(size(x,2),4);
        f     = {SBS,nSBS};
    case 'HBS'
        %hard block scaling
        HBS   = @(x) 1./nthroot(size(x,2),2);
        nHBS  = @(x) x./nthroot(size(x,2),2 );
        f     = {HBS,nHBS};
    case 'SHBS'
        %super hard block scaling
        SHBS  = @(x) 1./size(x,2);
        nSHBS = @(x) x./size(x,2);
        f     = {SHBS,nSHBS};
    case 'SBVS'
        %soft block variance scaling
        SBVS  = @(x) nthroot(size(x,2),4)./nthroot(sum(std(x,[],1).^2),2);
        nSBVS = @(x) x*nthroot(size(x,2),4)./nthroot(sum(std(x,[],1).^2),2);
        f     = {SBVS,nSBVS};
    case 'HBVS'
        %hard block variance scaling
        HBVS  = @(x) 1./power(sum(std(x,[],2).^2),0.5);
        nHBVS = @(x) x./power(sum(std(x,[],2).^2),0.5);
        f     = {HBVS,nHBVS};
    case 'SHBVS'
        %super hard block variance scaling
        SHBVS = @(x) 1./sum(std(x,[],2));
        nSHBVS= @(x) x./sum(std(x,[],2));
        f     = {SHBVS,nSHBVS};
    case 'SBRS'
        %soft block rank scaling
        SBVS  = @(x) nthroot(size(x,2),4)./nthroot(sum(std(x,[],1).^2),2);
        SBRS  = @(x) SBVS(x)*pseudoRank(x);
        nSBRS = @(x) x*SBVS(x)*pseudoRank(x);
        f     = {SBRS,nSBRS};
    case 'HBRS'
        %hard block rank scaling
        HBVS  = @(x) 1./power(sum(std(x,[],2).^2),0.5);
        HBRS  = @(x) HBVS(x)*pseudoRank(x);
        nHBRS = @(x) x*HBVS(x)*pseudoRank(x);
        f     = {HBRS,nHBRS};
    case 'SHBRS'
        %super hard block rank scaling
        SHBVS  = @(x) 1./sum(std(x,[],2));
        SHBRS  = @(x) SHBVS(x)*pseudoRank(x);
        nSHBRS = @(x) x*SHBVS(x)*pseudoRank(x);
        f      = {SHBRS,nSHBRS};
    otherwise
        keyboard
end
[trainSet,calSet,testSet]=applyBlockscaling(f,...
    {xPtrain0,xPtrain1,xPtrain2},{xStrain0,xStrain1,xStrain2},...
    {xPcalibration0,xPcalibration1,xPcalibration2},{xScalibration0,xScalibration1,xScalibration2},...
    {xPtest0,xPtest1,xPtest2},{xStest0,xStest1,xStest2});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HelperFunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trainSet,calSet,testSet]=applyBlockscaling(f,xPTr,xSTr,xPCa,xSCa,xPTe,xSTe)
f1=f{1};
f2=f{2};
xPtrain0=xPTr{1};
xPtrain1=xPTr{2};
xPtrain2=xPTr{3};
xStrain0=xSTr{1};
xStrain1=xSTr{2};
xStrain2=xSTr{3};
xPcalibration0=xPCa{1};
xPcalibration1=xPCa{2};
xPcalibration2=xPCa{3};
xScalibration0=xSCa{1};
xScalibration1=xSCa{2};
xScalibration2=xSCa{3};
xPtest0=xPTe{1};
xPtest1=xPTe{2};
xPtest2=xPTe{3};
xStest0=xSTe{1};
xStest1=xSTe{2};
xStest2=xSTe{3};
pTrainset=[{f2(xPtrain0)},{f2(xPtrain1)},{f2(xPtrain2)}];
sTrainset=[{f2(xStrain0)},{f2(xStrain1)},{f2(xStrain2)}];
pCalset=[{xPcalibration0.*f1(xPtrain0)},{xPcalibration1.*f1(xPtrain1)},{xPcalibration2.*f1(xPtrain2)}];
sCalset=[{xScalibration0.*f1(xStrain0)},{xScalibration1.*f1(xStrain1)},{xScalibration2.*f1(xStrain2)}];
pTestset=[{xPtest0.*f1(xPtrain0)},{xPtest1.*f1(xPtrain1)},{xPtest2.*f1(xPtrain2)}];
sTestset=[{xStest0.*f1(xStrain0)},{xStest1.*f1(xStrain1)},{xStest2.*f1(xStrain2)}];
%%
trainSet={pTrainset{1}, pTrainset{2}, pTrainset{3}, sTrainset{1}, sTrainset{2}, sTrainset{3}};
calSet={pCalset{1}, pCalset{2}, pCalset{3}, sCalset{1}, sCalset{2}, sCalset{3}};
testSet={pTestset{1}, pTestset{2}, pTestset{3}, sTestset{1}, sTestset{2}, sTestset{3}};
end