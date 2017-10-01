function [] = run_wtsi(inputFile, allOutputFile, outputDir, minNumberOfSignature, maxNumberOfSignature)

if ~exist('minNumberOfSignature', 'var')
    minNumberOfSignature = 3
end
if ~exist('maxNumberOfSignature', 'var')
    maxNumberOfSignature = 12
end

%% Open matlabpool
if ( matlabpool('size') == 0 )
    matlabpool open 12; % opens the default matlabpool, if it is not already opened
end

%% Define parameters
iterationsPerCore = 84;
stability = zeros(maxNumberOfSignature, 1);
reconstructionError = zeros(maxNumberOfSignature, 1);

%% Sequentially deciphering signatures between minNumberOfSignature and maxNumberOfSignature
for totalSignatures = minNumberOfSignature : maxNumberOfSignature
    
    % Decipher the signatures of mutational processes from catalogues of mutations
    [input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
        decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, ...
            [outputDir '/' num2str(totalSignatures) '_signatures.mat'] );
    
    % Record the stability and average Frobenius reconstruction error
    stability(totalSignatures-minNumberOfSignature+1) = mean(processStabAvg);
    reconstructionError(totalSignatures-minNumberOfSignature+1) = norm(input.originalGenomes - processes*exposures, 'fro');
    
end

save(allOutputFile);
