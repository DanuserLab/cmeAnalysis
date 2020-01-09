

function [DAS_all, DAS_Cluster] = dasAnalysis(data, pm, alt_dir, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('alt_dir', @(x) ischar(x));
ip.addParameter('A_max_all', 300, @isnumeric);
ip.addParameter('RandSampling', false, @islogical);
ip.addParameter('bin_A', 1, @isnumeric);
ip.addParameter('condition_name', 'Treated', @ischar);
ip.addParameter('Category', {'Ia'}, @iscell);
ip.addParameter('save_fig', false, @islogical);
ip.addParameter('calibrate', true, @islogical);
ip.addParameter('factor_given', [], @isnumeric);
ip.addParameter('calculate_factor', false, @islogical);
ip.addParameter('example', false, @islogical);
ip.addParameter('data_ratio', 1, @isnumeric);
ip.addParameter('Bandwidth', 0.05, @isnumeric);
ip.parse(data, pm, alt_dir, varargin{:});
%==================================================
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
A_max_all = ip.Results.A_max_all;
pm = ip.Results.pm;

bin_A = ip.Results.bin_A;
condition_name = ip.Results.condition_name;
Category = ip.Results.Category;
alt_dir = ip.Results.alt_dir;
data = ip.Results.data;
RandSampling = ip.Results.RandSampling;
num_condition = max(size(data));
n_ch = size(data{1}(1).channels,2);
calibrate = ip.Results.calibrate;
data_ratio = ip.Results.data_ratio;
calculate_factor = ip.Results.calculate_factor;
Bandwidth = ip.Results.Bandwidth;
factor_given = ip.Results.factor_given;

%==================================================
DAS_all = cell(num_condition,1);
Track_info_all = cell(num_condition,1);
master_dir = cell(num_condition,1);
fileExisting = cell(num_condition,1);
%==================================================
currentFolder = cd;

for i_condition = 1:num_condition
    %---------------------------------
    if isempty(alt_dir)
    if n_ch == 1
    master_dir{i_condition} = [data{i_condition}(1).source '..' filesep];
    else
    master_dir{i_condition} = [data{i_condition}(1).source '..' filesep '..' filesep];
    end
    else
        cd(alt_dir);
        if i_condition == 1
           mkdir('ctrl');
           master_dir{i_condition} = [alt_dir filesep 'ctrl' filesep];
        elseif i_condition == 2
           mkdir('treated');
           master_dir{i_condition} = [alt_dir filesep 'treated' filesep];  
        else
           mkdir(['treated',num2str(i_condition)]);
           master_dir{i_condition} = [alt_dir filesep ['treated',num2str(i_condition)] filesep];   
        end
    cd(master_dir{i_condition});  
    end
    %---------------------------------
    file_path = master_dir{i_condition};
    %---------------------------------
    if RandSampling == true
    fileExisting{i_condition} = (exist(fullfile(file_path, 'Track_info.mat'), 'file') == 2)...
              &(exist(fullfile(file_path, 'Track_info_samp.mat'), 'file') == 2);
    else
    fileExisting{i_condition} = (exist(fullfile(file_path, 'Track_info.mat'), 'file') == 2);
    end
    %---------------------------------
    if (fileExisting{i_condition} == false) || (pm.overwriteTrack_info == true)
    [Track_info] = dasRead(data{i_condition},master_dir{i_condition}, 'Category', Category);
    save([file_path filesep 'Track_info.mat'],'Track_info','-v7.3');
    Track_info_all{i_condition} = Track_info;
    clear Track_info;
    [Track_info_child] = dasReadChild(data{i_condition},master_dir{i_condition}, 'Category', Category);   
    save([file_path filesep 'Track_info_child.mat'],'Track_info_child','-v7.3');
    clear Track_info_child;
    else
    S = load([file_path filesep 'Track_info.mat']);
    Track_info_all{i_condition} = S.Track_info;
    clear S;
    disp('dasRead has been run ')
    end
end

%==========================================================================

if RandSampling == true
    if (fileExisting{i_condition} == false) || (pm.overwriteTrack_info == true)
        Track_info_all = dasRandSampling(data,master_dir,Track_info_all);
    else
        for i_condition = 1:num_condition
        cd(master_dir{i_condition});
        cd('Track_info');
        S = load('Track_info_samp.mat');
        Track_info_all{i_condition} = S.Track_info_samp;
        clear S;
        end
    end
end
%==========================================================================
if pm.i_max_D_adj == false
    A_max_all = pm.i_max_D;
else
    I_tem = [];
for i_condition = 1:num_condition
    for i_mov = 1:max(size(Track_info_all{i_condition}))  
        I_tem = [I_tem;max(Track_info_all{i_condition}{i_mov}.Max)];
    end   
end
    A_max_all = floor(max(I_tem)+0.5);
end

%==========================================================================
file_path = alt_dir;
DAS_done = (exist(fullfile(file_path, 'DAS_all.mat'), 'file') == 2);
%==========================================================================
if (pm.overwriteDAS == true) || (DAS_done == false)
%==========================================================================
    %---------------------------------
    if (DAS_done == false) || (calculate_factor == true)
    for i_condition = 1:num_condition
    disp('Processing raw data...')
    DAS_all{i_condition} = dasDataProcess(data{i_condition}, Track_info_all{i_condition});
    if i_condition == 1
        if calculate_factor == true
            DAS_all{1}.factor = dasCalibrate(data(1),DAS_all(1));
        else
            DAS_all{1}.factor = [1 0];
        end
    end
    end
    %---------------------------------
    else 
    disp('DAS calculation has been done before, reloading reusable data...')
        S = load([file_path filesep 'DAS_all.mat']);
        DAS_all_temp = S.DAS_all;
        clear S;
        factor_temp = DAS_all_temp{1}.factor;
        clear DAS_all_temp; 
    for i_condition = 1:num_condition
        disp('Processing raw data...')
        DAS_all{i_condition} = dasDataProcess(data{i_condition}, Track_info_all{i_condition});
        if i_condition == 1
            DAS_all{1}.factor = factor_temp;
        end
    end
    end
    %---------------------------------
for i_condition = 1:num_condition
    %---------------------------------
    if calibrate == true 
        factor_temp = DAS_all{1}.factor;
    else
        if isempty(factor_given)
        factor_temp = [1 0];
        else
            factor_temp = factor_given;
        end
    end
    if i_condition == 1
        [Dfunc,~, ~, ~] = dasComputeD(data{i_condition}, Track_info_all{i_condition}, [1 0], pm, A_max_all, 'bin_A', bin_A, 'data_ratio', data_ratio);
    elseif (i_condition > 1) && (pm.recalculateD == true)
        [Dfunc,~, ~, ~] = dasComputeD(data{i_condition}, Track_info_all{i_condition}, [1 0], pm, A_max_all, 'bin_A', bin_A, 'data_ratio', data_ratio);
    end

    DAS_all{i_condition} = dasComputeDas(data{i_condition}, Track_info_all{i_condition}, DAS_all{i_condition}, Dfunc, factor_temp, A_max_all, 'bin_A', bin_A);
    cd(master_dir{i_condition});
    %---------------------------------
end
%==========================================================================
save([file_path filesep 'DAS_all.mat'],'DAS_all','-v7.3');
%==========================================================================
else
    disp('DAS calculation has been finished')
    S = load([file_path filesep 'DAS_all.mat']);
    DAS_all = S.DAS_all;
    clear S;
%==========================================================================    
end
%==========================================================================

cd(currentFolder);

