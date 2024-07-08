function calibrate_classical_pls_manual(calibrationSet,yCalibration,randSeed)
if ~exist('randSeed','var')
    randSeed=1;
end
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
pCalset=calibrationSet{1};
sCalset=calibrationSet{4};
Q2Y=zeros(1,10);
R2Y=zeros(1,10);
xTemp=[pCalset,sCalset];
yTemp=yCalibration;
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,set1);
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    Calib=custom_pls(xTemp(set1,:),...
        yTemp(set1,:),10,...
        "isFirstblock",true,"needsPreprocess",false);
    Calib=nipalspredict(xTemp(set2,:),...
        Calib,'Y',yTemp(set2,:));
    Q2Y(j,:)=Calib.predict.Q2Y;
    R2Y(j,:)=Calib.Ry;
end
selectionHelperPlot(Q2Y,R2Y);
h=msgbox("Write down the results."+newline+...
    "close to continue");
while ishghandle(h)
    pause(1)
end
