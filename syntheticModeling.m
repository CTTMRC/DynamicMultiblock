calibrationTable=readtable('.\Data\metaParameter.csv');
%block scale to unit norm
SSN   = @(x) sqrt(sum(sum(x.^2)));
nSSN  = @(x) x./sqrt(sum(sum(x.^2)));
% ref listed in Data preprocessing for multiblock modelling –
% A systematization with new methods DOI:https://doi.org/10.1016/j.chemolab.2020.103959
SBS   = @(x) 1./nthroot(size(x,2),4);                                   %soft block scaling
HBS   = @(x) 1./nthroot(size(x,2),2);                                   %hard block scaling
SHBS  = @(x) 1./size(x,2);                                              %super hard block scaling
SVS   = @(x) nthroot(size(x,2),4)./nthroot(sum(std(x,[],1).^2),2);      %soft block variance scaling
HBVS  = @(x) 1./power(sum(std(x,[],2).^2),0.5);                         %hard block variance scaling
SHBVS = @(x) 1./sum(std(x,[],2));                                       %super hard block variance scaling
SBRS  = @(x) SVS(x)*pseudoRank(x);                                      %soft block rank scaling
HBRS  = @(x) HBVS(x)*pseudoRank(x);                                     %hard block rank scaling
SHBRS = @(x) SHBVS(x)*pseudoRank(x);                                    %super hard block rank scaling
nSBS  = @(x) x*1./nthroot(size(x,2),4);                              %soft block scaling
nHBS  = @(x) x*1./nthroot(size(x,2),2 );                              %hard block scaling
nSHBS = @(x) x*1./size(x,2);                                          %super hard block scaling
nSVS  = @(x) x*nthroot(size(x,2),4)./nthroot(sum(std(x,[],1).^2),2); %soft block variance scaling
nHBVS = @(x) x*1./power(sum(std(x,[],2).^2),0.5);                     %hard block variance scaling
nSHBVS= @(x) x*1./sum(std(x,[],2));                                   %super hard block variance scaling
nSBRS = @(x) x*SVS(x)*pseudoRank(x);                                  %soft block rank scaling
nHBRS = @(x) x*HBVS(x)*pseudoRank(x);                                 %hard block rank scaling
nSHBRS= @(x) x*SHBVS(x)*pseudoRank(x);                                %super hard block rank scaling
calibrationTable=rmmissing(calibrationTable);
rowNames=[repmat("C",[24,1]);repmat("S",[24,1])]+string(calibrationTable.R);
% rowNames(25)="-";
calibrationTable.Properties.RowNames=rowNames;
calibrationTable.R=[];
calibrationTable.Properties.VariableNames{4} = 'SOPLS-S';
calibrationTable.Properties.VariableNames{5} = 'SOPLS-P';
calibrationTable.Properties.VariableNames{6} = 'SODPLS-S';
calibrationTable.Properties.VariableNames{7} = 'SODPLS-P';
calibrationTable.Properties.VariableNames{8} = 'SOBDPLS-S1';
calibrationTable.Properties.VariableNames{9} = 'SOBDPLS-S2';
calibrationTable.Properties.VariableNames{10} = 'SOBDPLS-S3';
calibrationTable.Properties.VariableNames{11} = 'SOBDPLS-P1';
calibrationTable.Properties.VariableNames{12} = 'SOBDPLS-P2';
calibrationTable.Properties.VariableNames{13} = 'SOBDPLS-P3';
calibrationTable.Properties.VariableNames{14} = 'SODiPLS-S';
calibrationTable.Properties.VariableNames{15} = 'SODiPLS-P';
calibrationTableC=calibrationTable(1:24,:);
calibrationTableS=calibrationTable(25:48,:);
idx=[4,5,6,1,2,3];
dFF=fullfact([4,3,2]);
for i=1:size(dFF,1)
    weightCase      =dFF(i,1);
    blockWeight     =dFF(i,2);
    whichSpectra    =dFF(i,3);
    [trainSet,yTrain,~,~,testSet,yTest]=dataGen(weightCase,blockWeight,whichSpectra,'blockScaling','None');

    %%ClassicalPLS
    classicMdl=custom_pls(cat(2,[nSSN(trainSet{4}),nSSN(trainSet{1})]),...
        yTrain,calibrationTableC.Classic(i),...
        "isFirstblock",true,'needsPreprocess',true,'preprocessType','center');
    
    %%DPLS
    dplsMdl=custom_pls(cat(2,[nSSN(trainSet{4}),nSSN(trainSet{5}),nSSN(trainSet{6}),...
        nSSN(trainSet{1}),nSSN(trainSet{2}),nSSN(trainSet{3})])...
        ,yTrain,calibrationTableC.DPLS(i),...
        "isFirstblock",true,'needsPreprocess',true,'preprocessType','center');
    dplsMdl=getSuperScore(dplsMdl,{[1,141],[142,282],[283,423],[424,428],[429,433],[434,438]});
    %%DiPLS
    diplsMdl=DiPLS(cat(2,[nSSN(trainSet{4}),nSSN(trainSet{1})])...
        ,yTrain,calibrationTableC.DiPLS(i),calibrationTableS.DiPLS(i),...
        'needsPreprocess',true,'preprocessType','center','isAugmented',false);
    diplsMdl=getSuperScore(diplsMdl,{[1,141],[142,146]});
    %%SOPLS
    soplsMdl=SOPLS({trainSet{4},trainSet{1}}...
        ,yTrain,{calibrationTableC.("SOPLS-S")(i),calibrationTableC.("SOPLS-P")(i)}...
        ,'nipals','needsPreprocess',true,'preprocessType','center');
    
    %%SODPLS
    sodplsMdl=SOPLS({cat(2,[trainSet{4},trainSet{5},trainSet{6}]),...
        cat(2,[trainSet{1},trainSet{2},trainSet{3}])},...
        yTrain,{calibrationTableC.("SODPLS-S")(i),calibrationTableC.("SODPLS-P")(i)},...
        'nipals','needsPreprocess',true,'preprocessType','center');
    sodplsMdl(1)=getSuperScore(sodplsMdl(1),{[1,141],[142,282],[283,423]});
    sodplsMdl(2)=getSuperScore(sodplsMdl(2),{[1,5],[6,10],[11,15]});
    %%SOBDPLS
    steps=calibrationTableC{i,["SOBDPLS-S1","SOBDPLS-S2","SOBDPLS-S3","SOBDPLS-P1","SOBDPLS-P2","SOBDPLS-P3"]};
    ithStep=idx(steps~=0);
    nComp=steps(steps~=0);
    ithComponent=cell(size(ithStep));
    for j=1:length(nComp);ithComponent{j}=nComp(j); end
    sodaplsMdl=SOPLS(cat(2,trainSet(ithStep)),...
        yTrain,ithComponent...
        ,'nipals','needsPreprocess',true,'preprocessType','center');
    
    %%SODiPLS
    sodiplsMdl=SOPLS({trainSet{4},trainSet{1}},yTrain,...
        {calibrationTableC.("SODiPLS-S")(i),calibrationTableC.("SODiPLS-P")(i)},...
        'DiPLS',"dynamicOrder",[calibrationTableS.("SODiPLS-S")(i),calibrationTableS.("SODiPLS-P")(i)],...
        "isAugmented",false,'needsPreprocess',true,'preprocessType','center');
    
    %% Modeling - Testing
    classicTest=nipalspredict(cat(2,[testSet{4}/SSN(trainSet{4}),testSet{1}/SSN(trainSet{1})]),...
        classicMdl,'Y',yTest);
    dplsTest= nipalspredict(cat(2,[testSet{4}/SSN(trainSet{4}),testSet{5}/SSN(trainSet{5}),testSet{6}/SSN(trainSet{6}),...
        testSet{1}/SSN(trainSet{1}),testSet{2}/SSN(trainSet{2}),testSet{3}/SSN(trainSet{3})]),...
        dplsMdl,   'Y',yTest);
    diplsTest=DiPLSpredict(cat(2,[testSet{4}/SSN(trainSet{4}),testSet{1}/SSN(trainSet{1})]),...
        diplsMdl,  'Y',yTest);
    soplsTest=SOPLSpredict({testSet{4},testSet{1}},...
        soplsMdl,  'Y',yTest);
    sodplsTest=SOPLSpredict({cat(2,[testSet{4},testSet{5},testSet{6}]),...
        cat(2,[testSet{1},testSet{2},testSet{3}])},...
        sodplsMdl,  'Y',yTest);
    sodaplsTest=SOPLSpredict(cat(2,testSet(ithStep)),...
        sodaplsMdl,'Y',yTest);
    sodiplsTest=SOPLSpredict({testSet{4},testSet{1}},...
        sodiplsMdl,'Y',yTest);
    
    %%
    filename=".\Data\modelCollection_LVLs_"+string(dFF(i,1))+string(dFF(i,2))+string(dFF(i,3));
    save((filename),...
        "classicTest","classicMdl",...
        "dplsTest","dplsMdl",...
        "diplsTest","diplsMdl",...
        "soplsTest","soplsMdl",...
        "sodplsTest","sodplsMdl",...
        "sodaplsTest","sodaplsMdl",...
        "sodiplsTest","sodiplsMdl")
end