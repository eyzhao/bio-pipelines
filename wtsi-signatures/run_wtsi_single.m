function [] = run_wtsi_single(inputFile, outputFile, signatureCount)

%% Open matlabpool
if ( matlabpool('size') == 0 )
    matlabpool open 8; % opens the default matlabpool, if it is not already opened
end

%% Define parameters
iterationsPerCore = 125;

% Decipher the signatures of mutational processes from catalogues of mutations
decipherMutationalProcesses(iterationsPerCore, signatureCount, inputFile, outputFile);

