function coreData = runIterations( ...
    inputPath, outputPath, iterationsPerCore, numberProcessesToExtract ...
)

    %% Simple validation of the input params
   %%totalCores = matlabpool('size');
  
   if ( exist('iterationsPerCore', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of iterations that need to be performed for each lab!' );
   end
   
   if ( exist('numberProcessesToExtract', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of processes that need to be extracted!' );
   end
   
   if ( exist('inputPath', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify an input file!' );
   end
     
   if ( exist('outputPath', 'var') == 0 )
       error( 'decipherMutationalProcesses: Please specify an output file!' );
   end
   
   %% Load the input file and validate its fields
   data = load(inputPath);
   printProgress = 1;

   mutationTypesToRemoveSet = data.mutationTypesToRemoveSet;

   [W H genomeErrorsPar genomesReconstructedPar] = ...
      extract(data.genomes, iterationsPerCore, numberProcessesToExtract, printProgress);

   %% Saving to the output file
   save(outputPath, 'W', 'H', 'genomeErrorsPar', 'genomesReconstructedPar', 'mutationTypesToRemoveSet');
   
end
