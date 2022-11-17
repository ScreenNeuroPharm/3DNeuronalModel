% This script converts the peak train stored in ascii format (.txt) in
% Matlab sparse format (.mat). This is the correct file format to run the
% analysis using the developed scripts and Matlab functions.
% The same directories tree are mantained. Time stamps and peak 
% amplitudes are stored.
% 
%                 Paolo Massobrio - last update 7th May 2020
% Updated March 2022

clear

start_folder = uigetdir(pwd, 'Select the MAIN Peak Detection TXT folder');
out_folder = strcat(start_folder, '_mat');
mkdir(out_folder);
cd(start_folder);
prompt = {'Duration (s):','Frequency (Hz):'};
dlgtitle = 'Simulation Specifics';
dims = [1 35];
definput = {'60','10000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
simtime = str2double(answer{1,1});
freq = str2double(answer{2,1});
dur_sample = simtime * freq;
% d = dir; % main directory del peak train
% for j = 3:length(d)
%     subfold_name = d(j).name;
%     cd(subfold_name);
%     subfold_path = pwd;
    dd = dir; % main directory del peak train
    for k = 3:length(dd)
        filename = dd(k).name;
        % pt = load(filename);
        peak_train = zeros(dur_sample,1);
        pt = importdata(filename); 
        if ~isempty(pt)
            pt = reshape(pt, [size(pt,1) * size(pt, 2), 1])*freq;
            if sum(isnan(pt))>0
                pt(isnan(pt))=[];
            pt=int64(pt);
            peak_train(pt)=1; 
            end
        end
%         peak_train(pt(:,1))=pt(:,2);
        peak_train = sparse(peak_train);
        artifact = [];
        cd(out_folder);
%         if ~exist(subfold_name)
%             mkdir(subfold_name);
%         end
%         cd(subfold_name);
        filename_out = [filename(1:end),'.mat'];
        save(filename_out,'peak_train','artifact','-mat');
        cd(start_folder);
    end
    cd ..
% end
cd ..\..
clear all
disp('End Of Processing!');