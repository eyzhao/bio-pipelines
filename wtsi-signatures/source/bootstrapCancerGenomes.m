function bootstrapGenomes = bootstrapCancerGenomes( genomes )

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

   normGenomes = genomes ./ repmat( sum(genomes), size(genomes, 1), 1 );
   bootstrapGenomes = mnrnd( round(sum(genomes)'), normGenomes' )';
end