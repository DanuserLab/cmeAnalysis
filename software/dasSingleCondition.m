function [] = dasSingleCondition(pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('pm', @(x) isstruct(x));
ip.parse(pm, varargin{:});

pm = ip.Results.pm;

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
dasAnalysis(pm);
S_DAS_all = load([pm.dir_alt filesep 'DAS_all.mat']);
DAS_all = S_DAS_all.DAS_all;
clear S_DAS_all;
%%
dasStat(pm, DAS_all,'overwriteStat', true);
S_DAS_stat_idx = load([pm.dir_alt filesep 'DAS_stat_idx.mat']);
DAS_stat = S_DAS_stat_idx.DAS_stat;
idx = S_DAS_stat_idx.idx;
test_p = S_DAS_stat_idx.test_p;
clear S_DAS_stat_idx;
%%
if ~strcmp(pm.fig_disp_mod,'none') 
disp('reading tracks...');
S = cell(pm.num_condition,1);
for i_condition = 1:pm.num_condition
    if i_condition == 1
        file_name_temp = [pm.dir_alt filesep 'ctrl' filesep 'Track_info.mat'];
    elseif i_condition==2
        file_name_temp = [pm.dir_alt filesep 'treated' filesep 'Track_info.mat'];
    else
        file_name_temp = [pm.dir_alt filesep 'treated' num2str(i_condition) filesep 'Track_info.mat'];
    end
    S{i_condition} = load(file_name_temp);
end
disp('finished reading tracks');
clear i_condition;
clear file_name_temp;

if pm.plot_fig1 == true
figH = dasFigure1(DAS_all, S{1}.Track_info,DAS_stat,idx, pm,'i_condition',1);

figH = dasFigure1_3D(DAS_all,DAS_stat,idx, pm,figH,'i_condition',1);

mkdir(pm.fig_dir);
if pm.save_image == true
dasSavePlot(figH, pm.fig_dir, 1,'PaperPosition',pm.PaperPosition);
end

end

figH = dasPlotCondition(DAS_all, idx, DAS_stat,test_p,pm);

if pm.save_image == true
mkdir(pm.fig_dir);
dasSavePlot(figH, pm.fig_dir, 3, 'PaperPosition',pm.PaperPosition);
end
clear figH;
end
