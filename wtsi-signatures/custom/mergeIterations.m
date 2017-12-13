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

   iteration_one = load(paths{1});

   number_of_files = length(paths);
   totalIterationsPerCore = size(iteration_one.genomeErrorsPar, 3);
   numberProcessesToExtract = size(iteration_one.W, 2) / totalIterationsPerCore;
   totalGenomes = size(iteration_one.H, 2);
   totalMutationTypes = size(iteration_one.W, 1);

   Wall = zeros( totalMutationTypes, totalIterationsPerCore * numberProcessesToExtract * number_of_files );
   Hall = zeros( numberProcessesToExtract * totalIterationsPerCore * number_of_files, totalGenomes );
   genomeErrors = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * number_of_files );
   genomesReconstructed = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * number_of_files );

   disp('Computing Wall and Hall matrices')
   stepAll = numberProcessesToExtract * totalIterationsPerCore;
   step = 1;
   for startAll = 1 : stepAll : size(Wall, 2)
       disp(['Loading data from ' paths{step}])
       iteration_data = load(paths{step});
       endAll = startAll + stepAll - 1;
       Wall(:, startAll:endAll) = iteration_data.W;
       Hall(startAll:endAll, :) = iteration_data.H;
       step = step + 1;
   end

   disp('\nCalculating genome error matrices')
   step = 1;
   stepAll = totalIterationsPerCore;
   for startAll = 1 : stepAll : size( genomeErrors, 3 )
       disp(['Loading data from ' paths{step}])
       iteration_data = load(paths{step});
       endAll = startAll + stepAll - 1;
       genomeErrors(:, :, startAll:endAll) = iteration_data.genomeErrorsPar;
       genomesReconstructed(:, :, startAll:endAll) = iteration_data.genomesReconstructedPar;
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
    
   mutationTypesToRemoveSet = iteration_one.mutationTypesToRemoveSet;

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

    disp(['Wrote output to ' outputPath])
   
end
