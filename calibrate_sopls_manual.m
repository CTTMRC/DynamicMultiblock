function calibrate_sopls_manual(calibrationSet,yCalibration,randSeed)
if ~exist('randSeed','var')
    randSeed=1;
end
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
pCalset=calibrationSet{1};
sCalset=calibrationSet{4};
Q2Y=zeros(1,8);
R2Y=zeros(1,8);
yTemp=yCalibration;
for k=1%1:2
    trainSet=sCalset;%xTemp{k}
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,set1);
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    trainSet1={trainSet(set1,:)};
    trainSet2={trainSet(set2,:)};
    Calib=SOPLS(trainSet1,yTemp(set1,:),{8},'nipals','needsPreprocess',false);
    Calib=SOPLSpredict(trainSet2,Calib,'Y',yTemp(set2,:));
    Q2Y(j,:,k)=Calib(1).predict.Q2Y;
    R2Y(j,:,k)=Calib(1).Ry;
end
end

selectionHelperPlot(Q2Y,R2Y)
% uiwait(helpdlg('On the next window, pick which you want to plot.'))
% Have specify which plots he wants to make
calibration.nComp=[];
calibration.BlockOrder=[];
selectedIndexes = inputdlg({'nComp'},'Select Informations');%selectedIndexes = inputdlg({'nComp','BlockOrder'},'Select Informations');
if isempty(selectedIndexes)
    % User clicked Cancel.
    return;
end
calibration.nComp=str2double(selectedIndexes{1});
% calibration.BlockOrder=str2double(selectedIndexes{2});
blockOne=sCalset;%calibration.BlockOrder};
blockTwo=pCalset;%setdiff(1:2,calibration.BlockOrder)};
Q2Y=zeros(1,8);
R2Y=zeros(1,8);
RR2Y=zeros(1,8);
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,set1);
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    trainSet1={blockOne(set1,:),blockTwo(set1,:)};
    trainSet2={blockOne(set2,:),blockTwo(set2,:)};
    Calib=SOPLS(trainSet1,yTemp(set1,:),{calibration.nComp,8},'nipals','needsPreprocess',false);
    Calib=SOPLSpredict(trainSet2,Calib,'Y',yTemp(set2,:));
    Q2Y(j,:)=Calib(2).predict.Q2Y;
    R2Y(j,:)=Calib(2).Ry;
    RR2Y(j)=max(Calib(end).predict.R2Y);
end
selectionHelperPlot(Q2Y,R2Y)
h=helpdlg("Write down the results."+newline+ "R2= "+ num2str(mean(RR2Y))+newline+"(c)");
while ishghandle(h)
    pause(1)
end
end
