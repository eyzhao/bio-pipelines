function [] = runIterations( ...
    inputPath, outputPath, iterationsPerCore, numberProcessesToExtract ...
)

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

   iterationsPerCore = str2num(iterationsPerCore);
   numberProcessesToExtract = str2num(numberProcessesToExtract);

   %% Load the input file and validate its fields
   disp('loading data')
   data = load(inputPath);
   printProgress = 1;

   mutationTypesToRemoveSet = data.mutationTypesToRemoveSet;

   disp('starting extraction')

   [W H genomeErrorsPar genomesReconstructedPar] = ...
      extract(data.genomes, iterationsPerCore, numberProcessesToExtract, printProgress);

   %% Saving to the output file

   % OCTAVE save command
   save('-mat', outputPath, 'W', 'H', 'genomeErrorsPar', 'genomesReconstructedPar', 'mutationTypesToRemoveSet');

   % MATLAB save command
   %save(outputPath, 'W', 'H', 'genomeErrorsPar', 'genomesReconstructedPar', 'mutationTypesToRemoveSet');

   disp(['Wrote output to ' outputPath])
   
end
