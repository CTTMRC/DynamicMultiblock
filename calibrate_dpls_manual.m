function calibrate_dpls_manual(calibrationSet,yCalibration,randSeed)

if ~exist('randSeed','var')
    randSeed=1;
end

rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
Q2Y=zeros(length(permutationSet),25);
R2Y=zeros(length(permutationSet),25);
calSet=cat(2,calibrationSet{[4,5,6,1,2,3]});
for j=1:length(permutationSet)
    set1=[1,permutationSet(j,:)];
    set2=setdiff(1:10,permutationSet(j,:));
    set1=cat(1,blockSet{set1});
    set2=cat(1,blockSet{set2});
    Calib=custom_pls(calSet(set1,:)...
        ,yCalibration(set1,:),25,...
        "isFirstblock",true,"needsPreprocess",false);
    Calib=nipalspredict(calSet(set2,:),...
        Calib,'Y',yCalibration(set2,:));
    Q2Y(j,:)=Calib.predict.Q2Y;
    R2Y(j,:)=Calib.Ry;
    blockProgress(j,length(permutationSet),j==1)
end
selectionHelperPlot(Q2Y,R2Y);

h=msgbox("Write down the results."+newline+...
    "close to continue");
while ishghandle(h)
    pause(1)
end
