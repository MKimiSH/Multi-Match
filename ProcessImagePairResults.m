function ProcessImagePairResults
%%%%%%%%%%%%%%%%%%%%%%%%

baseDir = 'C:\outputs\';
runName ='Miko_problematics_5_iters';
runName ='Miko_problematics_5_iters_blur_2';
runName ='Miko_problematics_5_iters_blur_1';
runName ='Miko_for_paper_20_iters';
runName = 'Miko_for_paper_40_iters'; %  'ForRejection';
runName = 'Miko_for_paper_Gilad_3';
% runName ='Miko_problematics_5_iters_blur_1_only_boat';

outputDir = [baseDir runName '\'];

files = dir([outputDir '*.txt']);

datasetNames = {};

% gather by dataset
for i = 1 : length(files)
    name = files(i).name;
    datasetName = name(1:strfind(name,'__') - 1);
    res = strfind(datasetNames,datasetName);
    validInd = find(not(cellfun('isempty', res)));
    if (isempty(validInd))
        datasetNames{end+1} = datasetName;
        validInd = length(datasetNames);
    end
    'y';
end

numDatasets = length(datasetNames);
sadErrors = cell(1,numDatasets);
overlapErrors = cell(1,numDatasets);


% gather by dataset
for i = 1 : length(files)
    name = files(i).name;
    datasetName = name(1:strfind(name,'__') - 1);
    res = strfind(datasetNames,datasetName);
    validInd = find(not(cellfun('isempty', res)));
    
    secondImg = str2double(name(strfind(name,'img1_') + 8));
    iteration = str2double(name(strfind(name,'iter_') + 5 : strfind(name,'_TV1_') - 1));
    sadError = str2double(name(strfind(name,'SAD_A_') + 6 : strfind(name,'OVERLAP_A_') - 2));
    overlapError = str2double(name(strfind(name,'OVERLAP_A_') + 10 : strfind(name,'SAD_OPT_') - 2));
    
    sadErrors{validInd}(secondImg-1,iteration) = sadError;
    overlapErrors{validInd}(secondImg-1,iteration) = overlapError;
    
    'y';
end

% h1 = figure;
h2 = figure;
% set(h1,'name','SAD errors');
set(h2,'name',[runName ' - Overlap errors']);

rows = ceil(numDatasets/4);
numSamples = size(overlapErrors{1},2);


orderedDB = {'bark Z R','boat Z R','graffiti viewpoint','wall viewpoint',...
    'bikes blur', 'trees blur','ubc compression','light'};
prefixes = {'Zoom + Rotation (Bark)     ' , 'Zoom + Rotation (Boat)     ' , 'Viewpoint change (Graffiti)     ' , ...
    'Viewpoint change (Wall)     ' , 'Blur (Bikes)     ' , 'Blur (Trees)     ' , 'JPEG compression (UBC)     ' , 'Brightness change (Light)     '};
latexString = '';

for i = 1 : numDatasets
    datasetName = datasetNames{i};
    datasetName = strrep(datasetName,'_',' ');
    IndexC = strfind(orderedDB,datasetName);
    index = find(not(cellfun('isempty', IndexC)));
   
    %     Zoom + Rotation (Bark) 				& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%}	& \small{4\%}\\ \hline
    %     Zoom + Rotation (Boat)					& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     Viewpoint change (Graffiti)				& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     Viewpoint change (Wall)				& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     Blur (Bikes)							& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     Blur (Trees)							& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     JPEG compression (UBC)				& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    %     Brightness change (Light)				& \small{1\%} 	& \small{2\%} 	& \small{3\%} 	& \small{4\%} 	& \small{4\%}\\ \hline
    
    figure(h2);
    subplot (rows,4,index); hold on
    OVERLAPs = overlapErrors{i};
    
    percentages = 100*sum(OVERLAPs<=0.2,2)/size(OVERLAPs,2);
    latexString = sprintf('%s \\small{%s} & \\small{%.1f\\%%} 	&\\small{%.1f\\%%} 	& \\small{%.1f\\%%} 	& \\small{%.1f\\%%}	&\\small{%.1f\\%%}\\\\ \\hline \n',...
       latexString,  prefixes{index},percentages(1), percentages(2), percentages(3), percentages(4), percentages(5));
    latexString = strrep(latexString,'.0','');
  %  latexString = strrep(latexString,'\%','');
    
    level1 = sum(OVERLAPs<=0.1,2);
    level2 = sum((OVERLAPs>0.1)&(OVERLAPs<=0.2),2);
    level3 = sum((OVERLAPs>0.2)&(OVERLAPs<=0.5),2);
    level4 = sum(OVERLAPs>0.5,2);
    
    Y_overlap = [level1 level2 level3 level4 ];
    
    bar(Y_overlap,'stack')
    % L1 = legend('< 0.1','0.1-0.2','0.2-0.5','> 0.5');
    
    title(datasetName);
    axis([0 6 0 numSamples]);
    set(gca,'XTick', 1:5);
    xlabel('image pair')
    %         ylabel('number of samples');
    %         legend('SAD errors');
end

'x'











