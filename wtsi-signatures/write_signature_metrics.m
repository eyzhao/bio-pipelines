function [] = write_signature_metrics(outputDir, metricsFile)

files = dir([outputDir '/*.mat']);
numfiles = length(files);

stability = zeros(numfiles, 1);
reconstruction = zeros(numfiles, 1);
numSignatures = zeros(numfiles, 1);

for k = 1:numfiles
    disp(k)
    data = load([outputDir '/' files(k).name]);
    stability(k) = mean(data.processStabAvg); 
    reconstruction(k) = norm(data.input.originalGenomes - data.processes * data.exposures, 'fro');
    split = regexp(files(k).name, '_', 'split'); disp(split); numSignatures(k) = int32(str2double(split(1)));
end

dlmwrite(metricsFile, [numSignatures, stability, reconstruction], 'delimiter', '\t')
