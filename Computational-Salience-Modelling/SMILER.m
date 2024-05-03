% Code adapted by Maximilian Parker and Alex Muhl-Richardson from example code on the SMILER Github pages by Calden Wloka - comments by the original author are retained

%% Check if SMILER is installed

% This call checks to see if iSMILER has been added to the search path;
% this is used as a proxy to check and see if SMILER has already been set
% up for this particular MATLAB session.
% If you are creating your own scripts on a system for which iSMILER(true)
% has been run, you can omit this step.
if(exist('iSMILER.m', 'file') ~= 2)
    pathroot = mfilename('fullpath'); % get the current file location
    [pathroot, ~, ~] = fileparts(pathroot); % trim off the file name so we get the current directory
    cd('../../smiler_matlab_tools'); % navigate to where iSMILER is located
    iSMILER; % add SMILER models to the MATLAB path
    cd(pathroot); % return to our original location to execute this example
end

%% Set up the default experiment
models = {'LDS', 'GBVS', 'FES'}; % a cell array of the SMILER model codes which we wish to execute

% The next few lines create a dynamically allocated array of function
% handles to invoke the models specified in the previous line
modfun = cell(length(models),1);
for i = 1:length(models)
    modfun{i} = str2func([models{i}, '_wrap']);
end

input_set = dir('../input_images'); % get the list of images located in the example directory
input_set = input_set(3:end);  % trim folder navigation elements '.' and '..'

% if the output directory does not yet exist, make it
if(~exist('../output_maps_default', 'dir'))
    mkdir('../output_maps_default');
end

% check to see if the output directories exist to save the output maps for
% each model. If the directories do not exist, make them
for j = 1:length(models)
    if(~exist(['../output_maps_default/', models{j}], 'dir'))
        mkdir(['../output_maps_default/', models{j}]);
    end
end
    
%% Calculate default output maps and save them
% loop through images in the outer loop to save on imread commands
disp('Now starting the experiment using default parameters');
for i = 1:length(input_set)
    img = imread(['../input_images/', input_set(i).name]); % read in the image
    for j = 1:length(models)
        disp(['Executing model ', models{j}, ' on image ', num2str(i), ' of ', num2str(length(input_set))]);
        salmap = modfun{j}(img); % execute the jth model on the ith image
        imwrite(salmap, ['../output_maps_default/', models{j}, '/', input_set(i).name]); % save the saliency map
    end
end
disp(' '); % create a space in the display output before starting the next experiment

