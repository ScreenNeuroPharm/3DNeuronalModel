clear all
clc

%% Selection of the peak folder
start_folder = uigetdir(pwd, 'Select the PeakTrains_mat folder');
if isempty(strfind(start_folder, 'PeakTrains_mat'))
    f = errordlg('Folder not corrected', 'Folder Error');
    return
end
disp(start_folder);
cd(start_folder);
phase_dir = dir;

if length(phase_dir) == 2
    f = errordlg('There are not PeakTrain Files', 'Folder Error');
    return
end

cd ..

%% Prompt for data input
prompt = {'Number of cells:','Number of layers:'};
dlgtitle = 'Layer specifics';
dims = [1 35];
definput = {'1091', '6'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
n_cells = str2double(answer{1,1});
num_layer = str2double(answer{2,1});

%% Creation of output folder
cd(start_folder);

filename = strcat(start_folder,'_split');
if ~exist(filename)
    mkdir(filename);
end
cd(filename);
new_folder = pwd;

fn_l0_sx = strcat(new_folder,'\layer_0001_sx');
if ~exist(fn_l0_sx)
    mkdir(fn_l0_sx);
end
fn_l1_sx = strcat(new_folder,'\layer_0002_sx');
fn_l2_sx = strcat(new_folder,'\layer_0003_sx');
fn_l3_sx = strcat(new_folder,'\layer_0004_sx');
fn_l4_sx = strcat(new_folder,'\layer_0005_sx');
fn_l5_sx = strcat(new_folder,'\layer_0006_sx');

fn_l0_dx = strcat(new_folder,'\layer_0001_dx');
if ~exist(fn_l0_dx)
    mkdir(fn_l0_dx);
end
fn_l1_dx = strcat(new_folder,'\layer_0002_dx');
fn_l2_dx = strcat(new_folder,'\layer_0003_dx');
fn_l3_dx = strcat(new_folder,'\layer_0004_dx');
fn_l4_dx = strcat(new_folder,'\layer_0005_dx');
fn_l5_dx = strcat(new_folder,'\layer_0006_dx');

%% Copying peaktrain divided into layers
cd(start_folder);
% corrent = pwd;
peak_folder = dir;
for j = 3:length(peak_folder)
    el = split(peak_folder(j).name,'.');
    el = split(el(1),'_');
    pop = el{2};
    el = el{end};
    el_num = str2double(el);
    % cd(corrent); 
   
    if(strcmp(pop, 'dx'))
        fn_l0 = fn_l0_dx;
        fn_l1 = fn_l1_dx;
        fn_l2 = fn_l2_dx;
        fn_l3 = fn_l3_dx;
        fn_l4 = fn_l4_dx;
        fn_l5 = fn_l5_dx;
    else 
        fn_l0 = fn_l0_sx;
        fn_l1 = fn_l1_sx;
        fn_l2 = fn_l2_sx;
        fn_l3 = fn_l3_sx;
        fn_l4 = fn_l4_sx;
        fn_l5 = fn_l5_sx;
    end
    if el_num <= 1*n_cells
        copyfile(peak_folder(j).name,fn_l0);
    end
%     for ii = 2:num_layer
%         if el_num <= ii*n_cells && el_num > (ii-1)*n_cells && num_layer >= ii
%             if ~exist(fn_l1)      --> da gestiere nome
%                 mkdir(fn_l1);
%             end
%             copyfile(peak_folder(j).name,fn_l1);
%         end
%     end
    if el_num <= 2*n_cells && el_num > 1*n_cells && num_layer >= 2
        if ~exist(fn_l1)
            mkdir(fn_l1);
        end
        copyfile(peak_folder(j).name,fn_l1);
    end
    if el_num <= 3*n_cells && el_num > 2*n_cells && num_layer >= 3
        if ~exist(fn_l2)
            mkdir(fn_l2);
        end
        copyfile(peak_folder(j).name,fn_l2);
    end        
    if el_num <= 4*n_cells && el_num > 3*n_cells && num_layer >= 4
        if ~exist(fn_l3)
            mkdir(fn_l3);
        end
        copyfile(peak_folder(j).name,fn_l3);
    end
    if el_num <= 5*n_cells && el_num > 4*n_cells && num_layer >= 5
        if ~exist(fn_l4)
            mkdir(fn_l4);
        end
        copyfile(peak_folder(j).name,fn_l4);
    end
    if el_num <= 6*n_cells && el_num > 5*n_cells && num_layer >= 6
        if ~exist(fn_l5)
            mkdir(fn_l5);
        end
        copyfile(peak_folder(j).name,fn_l5);
    end
end
cd(start_folder);

EndOfProcessing (start_folder, 'Successfully accomplished');
clear
            