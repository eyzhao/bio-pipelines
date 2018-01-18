function [] = ...
    outputIterationData(inputPath, outputPath)
  
   %% Simple validation of the input params
   %%totalCores = matlabpool('size');
   if ( exist('inputPath', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify an input file!' );
   end
     
   if ( exist('outputPath', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify an input file!' );
   end
   
   %% Load the input file and validate its fields
   input = load(inputPath);
   if ( isfield(input, 'cancerType') == 0 || isfield(input, 'originalGenomes') == 0 ...
           || isfield(input, 'sampleNames') == 0 || isfield(input, 'subtypes') == 0 ...
           || isfield(input, 'types') == 0)
     error( 'decipherMutationalProcesses: Please specify an input file containing the variables: cancerType, originalGenomes, sampleNames, subtypes, types!' );
   end
   
   input.originalGenomes = double(input.originalGenomes)

   removeWeakMutationTypes = 0.01; % removes weak mutation types, i.e. reduces the dimmensions
   [genomes mutationTypesToRemoveSet] = removeWeak(input.originalGenomes, removeWeakMutationTypes);
     
   save(outputPath, 'genomes', 'mutationTypesToRemoveSet')
   disp(['File output to ' outputPath])
end
