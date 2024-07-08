function calibrate_sodapls_manual(calibrationSet,yCalibration,randSeed)
if ~exist('randSeed','var')
    randSeed=1;
end
rng(randSeed,'twister')
[blockSet,permutationSet]=cvBins(ones(size(calibrationSet{1})),10,'Reduced',true);
Q2Y=zeros(length(permutationSet),8,6);
R2Y=zeros(length(permutationSet),8,6);
xTemp=calibrationSet;
yTemp=yCalibration;
i=[4,5,6,1,2,3];
blockOne=xTemp(i(1));
blockTwo=[];
nComp={8};
iter=1;

for k=1:6 %First one is always going to be t0

    for j=1:length(permutationSet)
        set1=[1,permutationSet(j,:)];
        set2=setdiff(1:10,set1);
        set1=cat(1,blockSet{set1});
        set2=cat(1,blockSet{set2});
        trainSet1=cellfun(@(x)x(set1,:),[blockOne,blockTwo],'UniformOutput',false);
        trainSet2=cellfun(@(x)x(set2,:),[blockOne,blockTwo],'UniformOutput',false);
        
        Calib=SOPLS(trainSet1,yTemp(set1,:),nComp,'nipals','needsPreprocess',false);
        Calib=SOPLSpredict(trainSet2,Calib,'Y',yTemp(set2,:));
        Q2Y(j,:,k)=Calib(iter).predict.Q2Y;
        R2Y(j,:,k)=Calib(iter).Ry;
    end
    fig=selectionHelperPlot(Q2Y,R2Y);
    prompt = {'nComp','include?'};
    dlgtitle='Select Parameter to Continue';
    definput = {'0','1'};
    dims = [1 40];
    opts.Interpreter = 'tex';
    selectedIndexes = inputdlg(prompt,dlgtitle,dims,definput,opts);
    if isempty(selectedIndexes)
        % Clicked Cancel.
        includeFlag=0;
    else
        includeFlag=min(str2double(selectedIndexes{2}),str2double(selectedIndexes{1}));
    end
    if includeFlag&&k==1
        nComp{iter}=str2double(selectedIndexes{1});
        iter=iter+1;
    elseif includeFlag&&k>1
        blockOne=cat(2,[blockOne,xTemp(i(k))]);
        nComp{iter}=str2double(selectedIndexes{1});
        iter=iter+1;
    elseif~includeFlag
                
    else
        keyboard
    end
    nComp{iter}=8;
    if k<6
    blockTwo=xTemp(i(k+1));
    close(fig)
    end

end


setGlobalx(Q2Y)
h=helpdlg("Write down the results.");
while ishghandle(h)
    pause(1)
end
end