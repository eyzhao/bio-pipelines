function [Wall Hall genomeErrors genomesReconstructed idx idxS processes processesStd exposure exposureStd processStab processStabAvg clusterCompactness] = ... 
   deconvolute(genomesOriginal, totalIterationsPerCore, numberProcessesToExtract, processesDistance, removeWeakMutationTypes, totalCores)

%
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
%

   totalReplicates = 100;
   removeLastPercentage = 0.07;
   
   [genomes mutationTypesToRemoveSet] = removeWeak(genomesOriginal, removeWeakMutationTypes);
   totalMutationTypes = size(genomes, 1);
   totalGenomes = size(genomes, 2);
                                                      
   disp( ['Extracting ' num2str(numberProcessesToExtract) ' mutational signatures for ' num2str(totalIterationsPerCore) ' iterations on ' num2str(totalCores) ' labs:'] );
   
   spmd
       
    if ( labindex == 1 )
        printProgress = 1;
    else
        printProgress = 0;
    end
    
    [W H genomeErrorsPar genomesReconstructedPar] = ...
         extract(genomes, totalIterationsPerCore, numberProcessesToExtract, printProgress);
   end

   Wall = zeros( totalMutationTypes, numberProcessesToExtract * totalIterationsPerCore * length(W) );
   Hall = zeros( numberProcessesToExtract * totalIterationsPerCore * length(W), totalGenomes );
   genomeErrors = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * length(genomeErrorsPar) );
   genomesReconstructed = zeros(totalMutationTypes, totalGenomes, totalIterationsPerCore * length(genomesReconstructedPar) );

   stepAll = numberProcessesToExtract * totalIterationsPerCore;
   step = 1;
   for startAll = 1 : stepAll : size(Wall, 2)
       endAll = startAll + stepAll - 1;
       Wall(:, startAll:endAll) = W{step};
       Hall(startAll:endAll, :) = H{step};
       step = step + 1;
   end
   clear W H printProgress;

   step = 1;
   stepAll = totalIterationsPerCore;
   for startAll = 1 : stepAll : size( genomeErrors, 3 )
       endAll = startAll + stepAll - 1;
       genomeErrors(:, :, startAll:endAll) = genomeErrorsPar{step};
       genomesReconstructed(:, :, startAll:endAll) = genomesReconstructedPar{step};
       step = step + 1;
   end
   clear genomeErrorsPar genomesReconstructedPar averageIterationErrorPar;
   
   save( '/projects/ezhao_prj/temp/prelim_data.mat', ...
          'Wall', 'Hall', 'genomeErrors', 'genomes', ...
          'numberProcessesToExtract', ...
          'genomesReconstructed', 'removeLastPercentage');
   
    [ Wall Hall genomeErrors genomesReconstructed ] = ...
      filterOutIterations( Wall, Hall, genomeErrors, ...
                           numberProcessesToExtract, ...
                           genomesReconstructed, removeLastPercentage );


   [processes processesStd exposure exposureStd idx idxS processStab processStabAvg clusterCompactness] = ...
                                                          evaluateStability(Wall, Hall, ...
                                                          numberProcessesToExtract, ...
                                                          totalReplicates, ...
                                                          processesDistance);
    
   [processes processesStd Wall genomeErrors genomesReconstructed ] = ...
                                        addWeak( mutationTypesToRemoveSet, ... 
                                                 processes, processesStd, ...
                                                 Wall, genomeErrors, ...
                                                 genomesReconstructed );
    
end
