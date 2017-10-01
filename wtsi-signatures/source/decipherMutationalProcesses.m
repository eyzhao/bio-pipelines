function [input allProcesses allExposures idx processes exposures processStab processStabAvg] = decipherMutationalProcesses( totalIterationsPerCore, ...
                                                                                                                numberProcessesToExtract, ...
                                                                                                                 inputFileName, outputFileName)
% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.
   
   %% Define function specific constants
   processesDistance = 'cosine';
   removeWeakMutationTypes = 0.01; % removes weak mutation types, i.e. reduces the dimmensions

   %% Simple validation of the input params
   totalCores = matlabpool('size');
   if ( totalCores == 0)
       error( 'decipherMutationalProcesses: Please initialize a matlabpool!' );
   end
   
   if ( exist('totalIterationsPerCore', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of iterations that need to be performed for each lab!' );
   end
   
   if ( exist('numberProcessesToExtract', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of processes that need to be extracted!' );
   end
   
   if ( exist('inputFileName', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify an input file!' );
   end
     
   if ( exist('outputFileName', 'var') == 0 )
       error( 'decipherMutationalProcesses: Please specify an output file!' );
   end
   
   %% Load the input file and validate its fields
   input = load(inputFileName);
   if ( isfield(input, 'cancerType') == 0 || isfield(input, 'originalGenomes') == 0 ...
           || isfield(input, 'sampleNames') == 0 || isfield(input, 'subtypes') == 0 ...
           || isfield(input, 'types') == 0)
     error( 'decipherMutationalProcesses: Please specify an input file containing the variables: cancerType, originalGenomes, sampleNames, subtypes, types!' );
   end
   
   input.originalGenomes = double(input.originalGenomes)

   [allProcesses allExposures genomeErrors genomesReconstructed idx idxS processes processesStd exposures exposureStd processStab processStabAvg clusterCompactness] = ... ...
   deconvolute(input.originalGenomes, totalIterationsPerCore, numberProcessesToExtract, processesDistance, removeWeakMutationTypes, totalCores);
   
   %% Saving to the output file
   save(outputFileName, ...
        'allProcesses', 'allExposures', 'genomeErrors', 'genomesReconstructed', 'idx', ...
        'idxS', 'processes', 'processesStd', 'exposures', 'exposureStd', 'processStab', ...
        'processStabAvg', 'clusterCompactness', 'input');
   
end
