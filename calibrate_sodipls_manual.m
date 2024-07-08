function calibrate_sodipls_manual(calibrationSet,yCalibration,randSeed)
if ~exist('randSeed','var')
    randSeed=1;
end
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
Q2Y=zeros(length(permutationSet),8);
R2Y=zeros(length(permutationSet),8);
xTemp=calibrationSet([1,4]);
yTemp=yCalibration;
beta=nan(6,8,length(permutationSet),2);
blockOne=xTemp(2);
blockTwo=[];
nComp={8};
S=5;
for k=1:2
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,set1);
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    trainSet1=cellfun(@(x)x(set1,:),[blockOne,blockTwo],'UniformOutput',false);
    trainSet2=cellfun(@(x)x(set2,:),[blockOne,blockTwo],'UniformOutput',false);
    Calib=SOPLS(trainSet1,yTemp(set1,:),nComp,'DiPLS',"dynamicOrder",S,...
        "isAugmented",false,'needsPreprocess',false);
    Calib=SOPLSpredict(trainSet2,Calib,'Y',yTemp(set2,:));
    beta(:,1:size(Calib(k).T,2),j,k)=Calib(k).BetaX;
    R2Y(j,1:size(Calib(k).T,2),k)   =Calib(k).Ry;
    Q2Y(j,1:size(Calib(k).T,2),k)   =Calib(k).predict.Q2Y;
end
fig1=selectionHelperPlot(Q2Y,R2Y);
fig2=selectionHelperPlotBETA(beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibration.nComp=[];
calibration.dynS=[];
prompt = {'nComp','dynS'};
dlgtitle='Select Parameter to Continue';
definput = {'0','0'};
dims = [1 40];
opts.Interpreter = 'tex';
selectedIndexes = inputdlg(prompt,dlgtitle,dims,definput,opts);
if isempty(selectedIndexes)
    % Clicked Cancel.
    return
end
nComp{2}=8;
nComp{k}=str2double(selectedIndexes{1});
S(2)=5;
S(k)=str2double(selectedIndexes{2}); %#ok<AGROW>
blockTwo=xTemp{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if k==1
    close(fig1,fig2)
end
end
nComp=cell2mat(nComp);
setGlobalx(Q2Y)
setGlobaly(beta)

h=helpdlg("Write down the results."+newline+...
    "close to continue"+newline+"DO NOT CLICK OK"...
    +newline+"nComp= ["+num2str(nComp) + "]"...
    +newline+"S= ["+num2str(S) + "]");
while ishghandle(h)
    pause(1)
end
end