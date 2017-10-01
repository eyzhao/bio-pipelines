function [ Wall Hall genomeErrors genomesReconstructed ] = filterOutIterations( Wall, Hall, genomeErrors, ...
                                                                                numberProcessesToExtract, ...
                                                                                genomesReconstructed, removeLastPercentage )
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
    totalIterations = size(Wall, 2) / numberProcessesToExtract;
    totalRemoveIter = round(removeLastPercentage * totalIterations ); 
   
    closenessGenomes = zeros(totalIterations, 1);
    for i = 1 : totalIterations
        closenessGenomes(i) = norm(genomeErrors(:, :, i), 'fro');
    end
    
    [notUsed indexClosenessGenomes] = sort(closenessGenomes, 'descend');
    removeIterations = indexClosenessGenomes(1:totalRemoveIter);
   
    removeIterationSets = zeros( numberProcessesToExtract * totalRemoveIter, 1);
    
    for i = 1 : totalRemoveIter
        iStart = numberProcessesToExtract * ( removeIterations(i) - 1) + 1;
        iEnd = numberProcessesToExtract * removeIterations(i);
        removeIterationSets( (numberProcessesToExtract * ( i - 1) + 1):(numberProcessesToExtract * i) ) = ...
                             iStart:iEnd;
    end
   
    Wall( :, removeIterationSets ) = [];
    Hall( removeIterationSets, : ) = [];
    genomeErrors( :, :, removeIterations ) = [];
    genomesReconstructed( :, :, removeIterations ) = [];
    
end