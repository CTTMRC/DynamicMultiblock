function calibrate_sodpls_manual(calibrationSet,yCalibration,randSeed)
if ~exist('randSeed','var')
    randSeed=1;
end
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
pCalset=cat(2,calibrationSet{1:3});
sCalset=cat(2,calibrationSet{4:6});
Q2Y=zeros(length(permutationSet),8,2);
R2Y=zeros(length(permutationSet),8,2);
% xTemp= {pCalset,sCalset};
yTemp=yCalibration;
k=1;
trainSet=sCalset;%xTemp{k};
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


selectionHelperPlot(Q2Y,R2Y)
calibration.nComp=[];
calibration.BlockOrder=[];
selectedIndexes = inputdlg({'nComp'},'Select Informations');% inputdlg({'nComp','BlockOrder'},'Select Informations');
if isempty(selectedIndexes)
    % User clicked Cancel.
    return;
end
calibration.nComp=str2double(selectedIndexes{1});
% calibration.BlockOrder=str2double(selectedIndexes{2});
blockOne=sCalset;%xTemp{calibration.BlockOrder};
blockTwo=pCalset;%xTemp{setdiff(1:2,calibration.BlockOrder)};

RR2Y=zeros(1,8);
k=2;
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,set1);
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    trainSet1={blockOne(set1,:),blockTwo(set1,:)};
    trainSet2={blockOne(set2,:),blockTwo(set2,:)};
    Calib=SOPLS(trainSet1,yTemp(set1,:),{calibration.nComp,8},'nipals','needsPreprocess',false);
    Calib=SOPLSpredict(trainSet2,Calib,'Y',yTemp(set2,:));
    Q2Y(j,:,k)=Calib(2).predict.Q2Y;
    R2Y(j,:,k)=Calib(2).Ry;
    RR2Y(j)=max(Calib(end).predict.R2Y);
end
selectionHelperPlot(Q2Y,R2Y)
setGlobalx(Q2Y)
h=msgbox("Write down the results."+newline+...
    "close to continue");
while ishghandle(h)
    pause(1)
end
end

