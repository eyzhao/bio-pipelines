function [] = mergeIterations( ...
    originalInput, iterationPathsFile, outputPath ...
)

    %% Simple validation of the input params
   %%totalCores = matlabpool('size');
   
   if ( exist('originalInput', 'var') == 0 )
      error( 'mergeIterations: Please provide the file with original WTSI input data.' );
   end

   if ( exist('iterationPathsFile', 'var') == 0 )
      error( 'mergeIterations: Please provide a file with iteration file paths, one per line.' );
   end
     
   if ( exist('outputPath', 'var') == 0 )
       error( 'mergeIterations: Please specify an output file!' );
   end

   fid = fopen(iterationPathsFile);
   pathsData = textscan(fid, '%s', 'Delimiter', '\n');
   paths = pathsData{1};

   iterationData = cellfun(@(x) {load(x)}, paths);

   totalIterationsPerCore = size(iterationData{1}.genomeErrorsPar, 3);
   numberProcessesToExtract = size(iterationData{1}.W, 2) / totalIterationsPerCore;
   totalGenomes = size(iterationData{1}.H, 2);
   totalMutationTypes = size(iterationData{1}.W, 1);

   Wall = zeros( totalMutationTypes, totalIterationsPerCore * numberProcessesToExtract * length(iterationData) );
   Hall = zeros( numberProcessesToExtract * totalIterationsPerCore * length(iterationData), totalGenomes );
   genomeErrors = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * length(iterationData) );
   genomesReconstructed = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * length(iterationData) );

   stepAll = numberProcessesToExtract * totalIterationsPerCore;
   step = 1;
   for startAll = 1 : stepAll : size(Wall, 2)
       endAll = startAll + stepAll - 1;
       Wall(:, startAll:endAll) = iterationData{step}.W;
       Hall(startAll:endAll, :) = iterationData{step}.H;
       step = step + 1;
   end

   step = 1;
   stepAll = totalIterationsPerCore;
   for startAll = 1 : stepAll : size( genomeErrors, 3 )
       endAll = startAll + stepAll - 1;
       genomeErrors(:, :, startAll:endAll) = iterationData{step}.genomeErrorsPar;
       genomesReconstructed(:, :, startAll:endAll) = iterationData{step}.genomesReconstructedPar;
       step = step + 1;
   end

    % Required values from the main script are added right before they are needed here

   removeLastPercentage = 0.07;

    [ Wall Hall genomeErrors genomesReconstructed ] = ...
      filterOutIterations( Wall, Hall, genomeErrors, ...
                           numberProcessesToExtract, ...
                           genomesReconstructed, removeLastPercentage );

   processesDistance = 'cosine';
   totalReplicates = 100;

   [processes processesStd exposures exposureStd idx idxS processStab processStabAvg clusterCompactness] = ...
                                                          evaluateStability(Wall, Hall, ...
                                                          numberProcessesToExtract, ...
                                                          totalReplicates, ...
                                                          processesDistance);
    
   mutationTypesToRemoveSet = iterationData{1}.mutationTypesToRemoveSet;

   [processes processesStd Wall genomeErrors genomesReconstructed ] = ...
                                        addWeak( mutationTypesToRemoveSet, ... 
                                                 processes, processesStd, ...
                                                 Wall, genomeErrors, ...
                                                 genomesReconstructed );
 

    allProcesses = Wall;
    allExposures = Hall;
    input = load(originalInput);

   %% Saving to the output file
   save(outputPath, ...
        'allProcesses', 'allExposures', 'genomeErrors', 'genomesReconstructed', 'idx', ...
        'idxS', 'processes', 'processesStd', 'exposures', 'exposureStd', 'processStab', ...
        'processStabAvg', 'clusterCompactness', 'input');
   
end
