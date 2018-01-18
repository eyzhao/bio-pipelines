function [] = run_wtsi_single(inputFile, outputFile, signatureCount)

%% Open matlabpool

%% Define parameters
iterationsPerCore = 1000;

% Decipher the signatures of mutational processes from catalogues of mutations
decipherMutationalProcesses(iterationsPerCore, signatureCount, inputFile, outputFile, 1);

