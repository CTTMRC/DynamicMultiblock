function calibrate_dipls_manual(calibrationSet,yCalibration,randSeed)
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(yCalibration,10,'Reduced',true);
xTemp=[calibrationSet{1},calibrationSet{4}];
yTemp=yCalibration;
beta=zeros(7,8,length(permutationSet));
Q2Y=zeros(1,8);
R2Y=zeros(1,8);
for j=1:length(permutationSet)
    set1       =[1,permutationSet(j,:)];
    set2       =setdiff(1:10,permutationSet(j,:));
    set1       =cat(1,blockSet{set1});
    set2       =cat(1,blockSet{set2});
    Calib      =DiPLS(xTemp(set1,:),...
        yTemp(set1,:),8,6,...
        'needsPreprocess',false,'isAugmented',false);
    Calib=DiPLSpredict(xTemp(set2,:),...
        Calib,'Y',yTemp(set2,:));
    beta(:,:,j)=Calib.BetaX;
    Q2Y(j,:)   =Calib.predict.Q2Y;
    R2Y(j,:)   =Calib.Ry;
    blockProgress(j,length(permutationSet),j==1)
end
selectionHelperPlot(Q2Y,R2Y);
selectionHelperPlotBETA(beta);
h=helpdlg("Write down the results."+newline+...
    "close to continue"+newline+"DO NOT CLICK OK");
while ishghandle(h)
    pause(1)
end