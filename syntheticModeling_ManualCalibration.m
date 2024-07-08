dFF=fullfact([4,3,2]);
for i=1:size(dFF,1)
    %% Settings
    %%SEED
    randomSeed=1;
    randomGenerator='twister';
    %%RELATIVEWEIGHT
    weightCase=dFF(i,1);
    blockWeight=dFF(i,2);
    %%MISCELLANEA
    whichSpectra=dFF(i,3);
    noiseLevel=1;
    [~,~,standardCalset,yCalibration,~,~]=dataGen(weightCase,blockWeight,whichSpectra);
    calibrate_classical_pls_manual(standardCalset,yCalibration,1)
    calibrate_dpls_manual(standardCalset,yCalibration,1)
    calibrate_dipls_manual(standardCalset,yCalibration,1)
    calibrate_sopls_manual(standardCalset,yCalibration,1)
    calibrate_sodpls_manual(standardCalset,yCalibration,1)
    calibrate_sodapls_manual(standardCalset,yCalibration,1)
    calibrate_sodipls_manual(standardCalset,yCalibration,1)
end
