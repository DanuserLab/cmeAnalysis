if exist('data_all') == 0
    i_condition = 1;
    data_all{i_condition} = loadConditionData;
    i_condition = i_condition+1;
    add_cond = true;
    while add_cond == true
    x_temp = input('add more condition? (y/n)','s');
    if x_temp == 'y'
        data_all{i_condition} = loadConditionData;
        i_condition = i_condition+1;
    elseif x_temp == 'n'
        add_cond = false;
    else
        disp('input "y" or "n"');
    end
    end
end
clear i_condition;
clear x_temp;
clear add_cond;
num_condition = max(size(data_all));
%%
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
if exist('dir_alt') == 0
disp('Select an alternative folder to store DAS results:');
dir_alt = uigetdir(pwd, 'Select an alternative folder to store DAS results:');
end
%%
if exist('con_name') == 0
    con_name = cell(num_condition,1);
    for i_condition = 1:num_condition
        con_name{i_condition} = input(['Please input the name of condition # ', num2str(i_condition)],'s');
    end
end
%con_name{1} = 'WT';
clear i_condition;
%%
if exist('pm') == 0
    generate_pm = true;
    while generate_pm == true
    x_temp = input('use default parameters? (y/n)','s');
    if x_temp == 'y'
        pm = dasParameter(num_condition,dir_alt); 
        generate_pm = false;
    elseif x_temp == 'n'
        disp('please run "dasParameter" first');
        generate_pm = false;
    else
        disp('input "y" or "n"');
    end
    end
end
clear generate_pm;
clear x_temp;
%%
save([dir_alt filesep 'data_name_dir.mat'],'data_all','con_name','dir_alt')
%%
dasAnalysis(data_all,pm, dir_alt);
load([dir_alt filesep 'DAS_all.mat']);
%%
dasStat(dir_alt, pm, DAS_all,'overwriteStat', true);
load([dir_alt filesep 'DAS_stat_idx.mat']);
%%
disp('reading tracks...');
S = cell(num_condition,1);
for i_condition = 1:num_condition
    if i_condition == 1
        file_name_temp = [dir_alt filesep 'ctrl' filesep 'Track_info.mat'];
    elseif i_condition==2
        file_name_temp = [dir_alt filesep 'treated' filesep 'Track_info.mat'];
    else
        file_name_temp = [dir_alt filesep 'treated' num2str(i_condition) filesep 'Track_info.mat'];
    end
    S{i_condition} = load(file_name_temp);
end
disp('finished reading tracks');
clear i_condition;
clear file_name_temp;
%%
if pm.plot_fig1 == true
figH = dasFigure1(data_all, DAS_all, S{1}.Track_info,DAS_stat,idx, con_name,dir_alt,pm,'i_condition',1);

figH = dasFigure1_3D(data_all, DAS_all, S{1}.Track_info,DAS_stat,idx, Cx, con_name,dir_alt,pm,figH,'i_condition',1);

fig_dir = [dir_alt filesep 'figure'];
mkdir(fig_dir);
if pm.save_image == true
dasSavePlot(figH, pm.fig_dir, 1,'PaperPosition',pm.PaperPosition);
end

end
%%
figH = dasPlotCondition(data_all, DAS_all, idx, con_name,DAS_stat,test_p,pm,dir_alt);
%%
if pm.save_image == true
fig_dir = [dir_alt filesep 'figure'];
dasSavePlot(figH, pm.fig_dir, 3, 'PaperPosition',pm.PaperPosition);
end
clear figH;

