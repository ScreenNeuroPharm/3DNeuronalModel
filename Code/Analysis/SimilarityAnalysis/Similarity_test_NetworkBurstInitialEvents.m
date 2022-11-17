%%
clear
close all

%% Loading & folder management

[start_folder] = uigetdir(pwd,'Select the NetworkBurstDetectionFiles folder');
if isempty(strfind(start_folder, 'NetworkBurstDetectionFiles'))
    f = errordlg('Wrong folder', 'Folder Error');
    return
end
disp(start_folder)
cd(start_folder);
[StartFolderPath,StartFolderName] = fileparts(start_folder);
cd ../..
out_path = pwd; 
[ExpFolderPath,ExpFolderName] = fileparts(out_path);
out_folder = strcat(out_path, '\', ExpFolderName, '_SimilarityTest');
out_filename = strcat('SimilarityMatrix_',StartFolderName, '_dx.mat');
out_filename_fig = strcat('SimilarityMatrix_',StartFolderName, '_dx.fig');

if ~exist(out_folder)
    mkdir(out_folder);
end
cd(start_folder);

prompt = {'Duration (s):','PrepTime (s):','Frequency (Hz):'};
dlgtitle = 'Simulation Specifics';
dims = [1 35];
definput = {'60', '10', '10000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
simtime = str2double(answer{1,1});
freq = str2double(answer{3,1});
preptime = str2double(answer{2,1});

cd(start_folder);
d = dir;
initial_events = {};

for k = 3:length(d)
    % if isempty((strfind(d(k).name, 'parameters'))) && ~isempty((strfind(d(k).name, 'NetworkBurstDetection')))
    if isempty((strfind(d(k).name, 'parameters'))) && ~isempty((strfind(d(k).name, 'NetworkBurstDetection'))) && ~isempty((strfind(d(k).name, 'dx')))
        load(d(k).name);
        % initial_events = [initial_events, netBursts(:,1)- preptime*freq];
        initial_events = [initial_events, netBursts(:,1)];
    end
end


%% Computing VP distance

num = size(initial_events, 2);
dist = NaN(num);
Similarity_matrix = NaN(num);

for k = 1:size(initial_events,2)
    if ~isempty(initial_events{k})
        for j = k+1:size(initial_events,2)
            if ~isempty(initial_events{j})
                dist(k,j) = VP_compute_normalized_dist((initial_events{k}'/freq), (initial_events{j}'/freq), 0.5);
            else 
                dist(k,j) = NaN;
            end
        end
    else
        dist(k,k+1:end)= NaN;
    end
end
Similarity_matrix(~isnan(dist)) = 1-dist(~isnan(dist));

%% Safing and figures
        
cd(out_folder)
save(out_filename,'Similarity_matrix','-mat');

% Similarity_matrix(boolean(eye(size(Similarity_matrix))))= ones(size(Similarity_matrix, 2), 1);
heatmap(Similarity_matrix);
savefig(out_filename_fig);
cd(start_folder)

