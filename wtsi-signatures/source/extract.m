function [Wall Hall genomeErrors genomesReconstructed] = extract(genomes, totalIterations, totalProcesses, verbose)

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

    totalMutationTypes = size(genomes, 1);
    totalGenomes = size(genomes, 2);
    Wall = zeros( totalMutationTypes, totalProcesses * totalIterations );
    Hall = zeros( totalProcesses * totalIterations, totalGenomes );
    genomeErrors = zeros(totalMutationTypes, totalGenomes, totalIterations);
    genomesReconstructed = zeros(totalMutationTypes, totalGenomes, totalIterations);
    
    processCount = 1;
    for i = 1 : totalIterations
         
          % Generating boostrapped genomes
          bootstrapGenomes = max(bootstrapCancerGenomes( genomes ), eps);
          
          % Solving NMF for these boostrapped genomes
           [W, H] = nmf(bootstrapGenomes, totalProcesses, 0); % works good enough
          % Other NMF numerical methods
          % options = statset('TolX', 1e-12, 'TolFun', 1e-12, ...
          % 'MaxFunEvals', 100000000, 'MaxIter', 100000000);
          % [W, H] = nnmf(bootstrapGenomes, totalProcesses, 'algorithm', ...
          % 'mult', 'replicates', 1000); 
          % [W, H] = gdcls(bootstrapGenomes, totalProcesses, 100, 0.002, 'nonneg'); % sparse NMF
          % [W, H] = nmfsh_comb(bootstrapGenomes, totalProcesses, [10^-5 10], 0); % another good engine
 
        for j = 1 : totalProcesses
            total = sum( W(:, j) );
            W(:, j) = W(:, j) / total;
            H(j, :) = H(j, :) * total;
        end
        
        genomeErrors(:, :, i) = bootstrapGenomes -  W * H;
        genomesReconstructed(:, :, i) = W * H;
        Wall( :, processCount : (processCount + totalProcesses - 1) ) = W;
        Hall( processCount : (processCount + totalProcesses - 1), : ) = H;
        processCount = processCount + totalProcesses;
        
        if ( verbose == 1)
            disp( [ 'Iteration ' num2str(i) ' of ' num2str(totalIterations) ' for extracting ' num2str(totalProcesses) ...
                    ' mutational signatures in ' num2str(totalGenomes) ...
                    ' samples completed.'] );
        end
        
    end
    
end