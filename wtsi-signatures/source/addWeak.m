function [processes processesStd Wall genomeErrors genomesReconstructed ] = addWeak( mutationTypesToAddSet, ...
                                                                                                           processes_I, processesStd_I, ...
                                                                                                           Wall_I, genomeErrors_I, ...
                                                                                                           genomesReconstructed_I )
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

   totalMutTypes = size(Wall_I, 1) + length(mutationTypesToAddSet);
   processes = zeros( totalMutTypes, size(processes_I,2) );
   processesStd = zeros( totalMutTypes, size(processesStd_I,2) );
   Wall = zeros( totalMutTypes, size(Wall_I,2) );
   genomeErrors = zeros( totalMutTypes, size(genomeErrors_I, 2), size(genomeErrors_I, 3) );
   genomesReconstructed = zeros( totalMutTypes, size(genomesReconstructed_I, 2), size(genomesReconstructed_I, 3) );
   
   origArrayIndex = 1;
   for i = 1 : totalMutTypes
       
     if ( isempty( find(mutationTypesToAddSet == i) ) == 1 )
         
         processes(i, :) = processes_I(origArrayIndex, :);
         processesStd(i, :) = processesStd_I(origArrayIndex, :);
         Wall(i, :) = Wall_I(origArrayIndex, :);
         genomeErrors(i, :, :) = genomeErrors_I(origArrayIndex, :, :);
         genomesReconstructed(i, :, :) = genomesReconstructed_I(origArrayIndex, :, :);
         
         origArrayIndex = origArrayIndex + 1;
     end
     
   end    
end