% Define the directory and convert the structure to a table
dirPath = ".\Data\";
fileList = dir(dirPath);
fileTable = struct2table(fileList);

% Create a mask to filter the relevant files
fileMask = contains(fileTable.name, 'modelCollection_LVLs_') & ~contains(fileTable.name, 'TABLE');
fileTable = fileTable(fileMask, :);

% Extract directory and create the model list
firstDir = fileTable.folder{1, :};
modelList = strcat(firstDir, '\', fileTable.name);

% Initialize variables
step = 0;
factor1 = ["[0,0]"; "[1,1]"; "[0,2]"; "[2,0]"];
factor2 = ["[0.5,0.5]"; "[0.8,0.2]"; "[0.2,0.8]"];
factor3 = ["[1]"; "[2345]"];
span = 0:3:72;
fullBlocks = [];
fulDynOrd = [];
dynCoeff = 0.5;

% Initialize result tables
numRows = 24;
numVars = 7;
varNames = {'PLS', 'DPLS', 'DiPLS', 'SOPLS', 'SODPLS', 'SOBDPLS', 'SODiPLS'};
fullResultsR2 = table('Size', [numRows, numVars], 'VariableTypes', repmat("double", [1, numVars]), ...
    'VariableNames', varNames, 'RowNames', string(1:numRows)' + repmat("_R2  ", numRows, 1));
fullResultsRMSE = table('Size', [numRows, numVars], 'VariableTypes', repmat("double", [1, numVars]), ...
    'VariableNames', varNames, 'RowNames', string(1:numRows)' + repmat("_RMSE", numRows, 1));
fullResultsMAE = table('Size', [numRows, numVars], 'VariableTypes', repmat("double", [1, numVars]), ...
    'VariableNames', varNames, 'RowNames', string(1:numRows)' + repmat("_MAE ", numRows, 1));

% Loop through each model file and process
for i = 1:length(modelList)
    modelFile = modelList{i};
    load(modelFile);
    step = step + 1;

    % Initialize the temporary table for storing metrics
    tt = table('Size', [7, 3], 'VariableTypes', repmat("double", [1, 3]), 'VariableNames', {'R2_P'; 'RMSE_P'; 'MAE_P'}, ...
        'RowNames', {'PLS', 'DPLS', 'DiPLS', 'SOPLS', 'SODPLS', 'SODAPLS', 'SODiPLS'});

    % Extract metrics from different test results and store in the temporary table
    tt{'PLS', :} = [classicTest(end).predict.R2Y, classicTest(end).predict.RMSE, classicTest(end).predict.MAE];
    tt{'DPLS', :} = [dplsTest(end).predict.R2Y, dplsTest(end).predict.RMSE, dplsTest(end).predict.MAE];
    tt{'DiPLS', :} = [diplsTest(end).predict.R2Y, diplsTest(end).predict.RMSE, diplsTest(end).predict.MAE];
    tt{'SOPLS', :} = [soplsTest(end).predict.R2Y, soplsTest(end).predict.RMSE, soplsTest(end).predict.MAE];
    tt{'SODPLS', :} = [sodplsTest(end).predict.R2Y, sodplsTest(end).predict.RMSE, sodplsTest(end).predict.MAE];
    tt{'SODAPLS', :} = [sodaplsTest(end).predict.R2Y, sodaplsTest(end).predict.RMSE, sodaplsTest(end).predict.MAE];
    tt{'SODiPLS', :} = [sodiplsTest(end).predict.R2Y, sodiplsTest(end).predict.RMSE, sodiplsTest(end).predict.MAE];

    % Store the results in the main result tables
    fullResultsR2{step, :} = tt.R2_P';
    fullResultsRMSE{step, :} = tt.RMSE_P';
    fullResultsMAE{step, :} = tt.MAE_P';

    % Extract the model name and update row names
    nameExtraction = char(extract(modelFile, digitsPattern(3)));
    factor1Selector = str2double(nameExtraction(1));
    factor2Selector = str2double(nameExtraction(2));
    factor3Selector = str2double(nameExtraction(3));
    name = factor1(factor1Selector) + "-" + factor2(factor2Selector) + "-" + factor3(factor3Selector);
    fullResultsR2.Properties.RowNames{step} = cellstr(name + "_R2  ");
    fullResultsRMSE.Properties.RowNames{step} = cellstr(name + "_RMSE");
    fullResultsMAE.Properties.RowNames{step} = cellstr(name + "_MAE ");

    % Additional calculations for multi-block statistics
    blocks = calculateBlocks(dplsTest, diplsTest, soplsTest, sodplsTest, sodaplsTest, sodiplsTest, name);
    fullBlocks = [fullBlocks; blocks];
    
    dynOrder = calculateDynamicOrder(dplsTest, diplsTest, sodplsTest, sodaplsTest, sodiplsTest, name, dynCoeff);
    fulDynOrd = [fulDynOrd; dynOrder];
end

% Save the results
save(".\Data\TabularResults.mat", "fullResultsR2", "fullResultsRMSE", "fullResultsMAE");

% Calculate and save block results
blockDistance = calculateBlockDistance(fullBlocks);
blockMean = calculateBlockMean(fullBlocks);
blockSTD = calculateBlockSTD(fullBlocks);

% Calculate and save dynamic order results
dynDistance = calculateDynDistance(fulDynOrd);
dynMean = calculateDynMean(fulDynOrd);
dynSTD = calculateDynSTD(fulDynOrd);

% Save block and dynamic order results to Excel
% writetable(blockResults, "./Data/blockResults.xlsx", 'WriteRowNames', true, 'WriteVariableNames', true);
% writetable(dynResults, "./Data/dynResults.xlsx", 'WriteRowNames', true, 'WriteVariableNames', true);


%% Helper Functions

function blocks = calculateBlocks(dplsTest, diplsTest, soplsTest, sodplsTest, sodaplsTest, sodiplsTest, name)
    % Helper function to calculate blocks
    sodaplsSblock=0;
    sodaplsPblock=0;
    for j=1:size(sodaplsMdl,2)
        blockCheck=size(sodaplsTest(j).W,1);
        if blockCheck==141
            sodaplsSblock=sodaplsSblock+(sodaplsTest(j).predict.R2Y./sodaplsTest(end).predict.R2Y);
        elseif blockCheck==5
            sodaplsPblock=sodaplsPblock+(sodaplsTest(j).predict.R2Y./sodaplsTest(end).predict.R2Y);
        else
            keyboard
        end
        
    end
    blocks = table('Size', [2, 6], 'VariableTypes', repmat("double", [1, 6]), 'VariableNames', ...
        {'DPLS', 'DiPLS', 'SOPLS', 'SODPLS', 'SODAPLS', 'SODiPLS'}, 'RowNames', {name+"_S", name+"_P"});
    blocks{:,'DPLS'} = [sum(dplsTest.multiblockStats.tbRy(1:3)); sum(dplsTest.multiblockStats.tbRy(4:6))] / ...
        sum(dplsTest.multiblockStats.tbRy);
    blocks{:,'DiPLS'} = diplsTest.multiblockStats.tbRy / sum(diplsTest.multiblockStats.tbRy);
    blocks{1,'SOPLS'} = soplsTest(1).predict.R2Y / soplsTest(end).predict.R2Y;
    blocks{2,'SOPLS'} = soplsTest(2).predict.R2Y / soplsTest(end).predict.R2Y;
    blocks{1,'SODPLS'} = sodplsTest(1).predict.R2Y / sodplsTest(end).predict.R2Y;
    blocks{2,'SODPLS'} = sodplsTest(2).predict.R2Y / sodplsTest(end).predict.R2Y;
    blocks{1,'SODiPLS'} = sodiplsTest(1).predict.R2Y / sodiplsTest(end).predict.R2Y;
    blocks{2,'SODiPLS'} = sodiplsTest(2).predict.R2Y / sodiplsTest(end).predict.R2Y;
    blocks{1,'SODAPLS'} = sodaplsSblock;
    blocks{2,'SODAPLS'} = sodaplsPblock;
end

function dynOrder = calculateDynamicOrder(dplsTest, diplsTest, sodplsTest, sodaplsTest, sodiplsTest, name, dynCoeff)
    % Helper function to calculate dynamic order
    sodaplsSdynOrd=0;
    sodaplsPdynOrd=0;
    for j=1:size(sodaplsMdl,2)
        blockCheck=size(sodaplsTest(j).W,1);
        if blockCheck==141
            sodaplsSdynOrd=sodaplsSdynOrd+1;
        elseif blockCheck==5
            sodaplsPdynOrd=sodaplsPdynOrd+1;
        else
            keyboard
        end
        
    end
    
    
    dynOrder = table('Size', [2, 5], 'VariableTypes', repmat("double", [1, 5]), 'VariableNames', ...
        {'DPLS', 'DiPLS', 'SODPLS', 'SODAPLS', 'SODiPLS'}, 'RowNames', {name+"_S", name+"_P"});
    dynOrder{1,'SODAPLS'} = sodaplsSdynOrd - 1;
    dynOrder{2,'SODAPLS'} = sodaplsPdynOrd - 1;
    dynOrder{1,'SODPLS'}=sum(sodplsTest(1).multiblockStats.tbRy./max(sodplsTest(1).multiblockStats.tbRy)>dynCoeff);
    dynOrder{2,'SODPLS'}=sum(sodplsTest(2).multiblockStats.tbRy./max(sodplsTest(2).multiblockStats.tbRy)>dynCoeff);
    dynOrder{1,'SODiPLS'}=size(sodiplsTest(1).T_s,3)-1;
    dynOrder{2,'SODiPLS'}=size(sodiplsTest(2).T_s,3)-1;
    dynOrder{1,'DiPLS'}=size(diplsTest.T_s,3)-1;
    dynOrder{2,'DiPLS'}=size(diplsTest.T_s,3)-1;
    dynOrder{1,'DPLS'}=sum(dplsTest.multiblockStats.tbRy(1:3)./max(dplsTest.multiblockStats.tbRy(1:3))>dynCoeff);
    dynOrder{2,'DPLS'}=sum(dplsTest.multiblockStats.tbRy(4:6)./max(dplsTest.multiblockStats.tbRy(4:6))>dynCoeff);
    
end

function blockDistance=calculateBlockDistance(fullBlocks)
blocks55=contains(fullBlocks.Properties.RowNames,"[0.5,0.5]");
blocks82=contains(fullBlocks.Properties.RowNames,"[0.8,0.2]");
blocks28=contains(fullBlocks.Properties.RowNames,"[0.2,0.8]");
PP=contains(fullBlocks.Properties.RowNames,"_P");
SS=contains(fullBlocks.Properties.RowNames,"_S");
blockDistance=table('Size',[3,6],'VariableTypes',repmat("doublenan",[1,6]),'VariableNames',...
    fullBlocks.Properties.VariableNames,...
    'RowNames' ,cellstr(["[0.5,0.5]";"[0.8,0.2]";"[0.2,0.8]"]));
for i=1:6
    temp55=[fullBlocks{blocks55&SS,i},fullBlocks{blocks55&PP,i}];
        blockDistance{1,i}=mean(pdist2(temp55,[0.5,0.5],'euclidean'));

    temp82=[fullBlocks{blocks82&SS,i},fullBlocks{blocks82&PP,i}];
        blockDistance{2,i}=mean(pdist2(temp82,[0.8,0.2],'euclidean'));

    temp28=[fullBlocks{blocks28&SS,i},fullBlocks{blocks28&PP,i}];
        blockDistance{3,i}=mean(pdist2(temp28,[0.2,0.8],'euclidean'));    
end

end

function blockSTD=calculateBlockSTD(fullBlocks)
blocks55=contains(fullBlocks.Properties.RowNames,"[0.5,0.5]");
blocks82=contains(fullBlocks.Properties.RowNames,"[0.8,0.2]");
blocks28=contains(fullBlocks.Properties.RowNames,"[0.2,0.8]");
PP=contains(fullBlocks.Properties.RowNames,"_P");
SS=contains(fullBlocks.Properties.RowNames,"_S");
blockSTD=table('Size',[3,6],'VariableTypes',repmat("doublenan",[1,6]),'VariableNames',...
    fullBlocks.Properties.VariableNames,...
    'RowNames' ,cellstr(["[0.5,0.5]";"[0.8,0.2]";"[0.2,0.8]"]));
for i=1:6
    temp55=[fullBlocks{blocks55&SS,i},fullBlocks{blocks55&PP,i}];
    blockSTD{1,i}=std(temp55(:,1));
        
    temp82=[fullBlocks{blocks82&SS,i},fullBlocks{blocks82&PP,i}];
    blockSTD{2,i}=std(temp82(:,1));
       
    temp28=[fullBlocks{blocks28&SS,i},fullBlocks{blocks28&PP,i}];
    blockSTD{3,i}=std(temp28(:,1));
end

end

function blockMean=calculateBlockMean(fullBlocks)
blocks55=contains(fullBlocks.Properties.RowNames,"[0.5,0.5]");
blocks82=contains(fullBlocks.Properties.RowNames,"[0.8,0.2]");
blocks28=contains(fullBlocks.Properties.RowNames,"[0.2,0.8]");
PP=contains(fullBlocks.Properties.RowNames,"_P");
SS=contains(fullBlocks.Properties.RowNames,"_S");
blockMean=table('Size',[3,12],'VariableTypes',repmat("doublenan",[1,12]),'VariableNames',...
    [fullBlocks.Properties.VariableNames+"_S",fullBlocks.Properties.VariableNames+"_P"],...
    'RowNames' ,cellstr(["[0.5,0.5]";"[0.8,0.2]";"[0.2,0.8]"]));
for i=1:6
    temp55=[fullBlocks{blocks55&SS,i},fullBlocks{blocks55&PP,i}];
    blockMean{1,i}=mean(temp55(:,1));
    blockMean{1,i+6}=mean(temp55(:,2));

    temp82=[fullBlocks{blocks82&SS,i},fullBlocks{blocks82&PP,i}];
    blockMean{2,i}=mean(temp82(:,1));
    blockMean{2,i+6}=mean(temp82(:,2));

    temp28=[fullBlocks{blocks28&SS,i},fullBlocks{blocks28&PP,i}];
    blockMean{3,i}=mean(temp28(:,1));
    blockMean{3,i+6}=mean(temp28(:,2));
    
end
end

function dynDistance=calculateDynDistance(fulDynOrd)
dyn00=contains(fulDynOrd.Properties.RowNames,"[0,0]");
dyn11=contains(fulDynOrd.Properties.RowNames,"[1,1]");
dyn20=contains(fulDynOrd.Properties.RowNames,"[2,0]");
dyn02=contains(fulDynOrd.Properties.RowNames,"[0,2]");
PP=contains(fulDynOrd.Properties.RowNames,"_P");
SS=contains(fulDynOrd.Properties.RowNames,"_S");
dynDistance=table('Size',[4,5],'VariableTypes',repmat("doublenan",[1,5]),'VariableNames',...
    fulDynOrd.Properties.VariableNames,...
    'RowNames' ,cellstr(["[0,0]";"[1,1]";"[2,0]";"[0,2]"]));
for i=1:5%fulDynOrd.Properties.VariableNames
    temp00=[fulDynOrd{dyn00&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn00&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynDistance{1,i}=mean(pdist2(temp00,[0,0],'euclidean'));
    
    temp11=[fulDynOrd{dyn11&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn11&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynDistance{2,i}=mean(pdist2(temp11,[1,1],'euclidean'));
    
    temp20=[fulDynOrd{dyn20&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn20&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynDistance{3,i}=mean(pdist2(temp20,[2,0],'euclidean'));
    
    temp02=[fulDynOrd{dyn02&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn02&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynDistance{4,i}=mean(pdist2(temp02,[0,2],'euclidean'));
    
end
end

function dynMean=calculateDynMean(fulDynOrd)
dyn00=contains(fulDynOrd.Properties.RowNames,"[0,0]");
dyn11=contains(fulDynOrd.Properties.RowNames,"[1,1]");
dyn20=contains(fulDynOrd.Properties.RowNames,"[2,0]");
dyn02=contains(fulDynOrd.Properties.RowNames,"[0,2]");
PP=contains(fulDynOrd.Properties.RowNames,"_P");
SS=contains(fulDynOrd.Properties.RowNames,"_S");
dynMean=table('Size',[4,10],'VariableTypes',repmat("doublenan",[1,10]),'VariableNames',...
    [fulDynOrd.Properties.VariableNames+"_S",fulDynOrd.Properties.VariableNames+"_P"],...
    'RowNames' ,cellstr(["[0,0]";"[1,1]";"[2,0]";"[0,2]"]));
for i=1:5%fulDynOrd.Properties.VariableNames
    temp00=[fulDynOrd{dyn00&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn00&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynMean{1,i}=(mean((temp00(:,1))));
    dynMean{1,i+5}=(mean((temp00(:,2))));
    
    temp11=[fulDynOrd{dyn11&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn11&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynMean{2,i}=mean((temp11(:,1)));
    dynMean{2,i+5}=mean((temp11(:,2)));
    
    temp20=[fulDynOrd{dyn20&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn20&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynMean{3,i}=mean((temp20(:,1)));
    dynMean{3,i+5}=mean((temp20(:,2)));
    
    temp02=[fulDynOrd{dyn02&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn02&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynMean{4,i}=mean((temp02(:,1)));
    dynMean{4,i+5}=mean((temp02(:,2)));    
    
end
end

function dynSTD=calculateDynSTD(fulDynOrd)
dyn00=contains(fulDynOrd.Properties.RowNames,"[0,0]");
dyn11=contains(fulDynOrd.Properties.RowNames,"[1,1]");
dyn20=contains(fulDynOrd.Properties.RowNames,"[2,0]");
dyn02=contains(fulDynOrd.Properties.RowNames,"[0,2]");
PP=contains(fulDynOrd.Properties.RowNames,"_P");
SS=contains(fulDynOrd.Properties.RowNames,"_S");
dynSTD=table('Size',[4,10],'VariableTypes',repmat("doublenan",[1,10]),'VariableNames',...
    [fulDynOrd.Properties.VariableNames+"_S",fulDynOrd.Properties.VariableNames+"_P"],...
    'RowNames' ,cellstr(["[0,0]";"[1,1]";"[2,0]";"[0,2]"]));
for i=1:5%fulDynOrd.Properties.VariableNames
    temp00=[fulDynOrd{dyn00&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn00&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynSTD{1,i}=(std((temp00(:,1))));
    dynSTD{1,i+5}=(std((temp00(:,2))));
    
    temp11=[fulDynOrd{dyn11&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn11&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynSTD{2,i}=std((temp11(:,1)));
    dynSTD{2,i+5}=std((temp11(:,2)));
    
    temp20=[fulDynOrd{dyn20&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn20&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynSTD{3,i}=std((temp20(:,1)));
    dynSTD{3,i+5}=std((temp20(:,2)));
    
    temp02=[fulDynOrd{dyn02&SS,fulDynOrd.Properties.VariableNames(i)},fulDynOrd{dyn02&PP,fulDynOrd.Properties.VariableNames(i)}];
    dynSTD{4,i}=std((temp02(:,1)));
    dynSTD{4,i+5}=std((temp02(:,2)));    
    
end
end