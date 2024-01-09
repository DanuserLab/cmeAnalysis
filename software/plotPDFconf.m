function [fig_pdf] = plotPDFconf(data1, data2, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data1', @(x) iscell(x));
ip.addRequired('data2', @(x) iscell(x));
ip.addParameter('BoundedSupport', [4, Inf], @isnumeric); % BoundedSupport = [] (default) or [4, Inf] or [0, Inf]
ip.addParameter('save', false, @islogical);
ip.addParameter('plot1', false, @islogical);
ip.addParameter('varName', 'x', @ischar);
ip.addParameter('CondName', {'cond1','cond2'}, @iscell);
ip.addParameter('num_condition', 2, @isnumeric);
ip.parse(data1, data2, varargin{:});

save = ip.Results.save;
data1 = ip.Results.data1;
data2 = ip.Results.data2;
BoundedSupport = ip.Results.BoundedSupport;
varName = ip.Results.varName;
CondName = ip.Results.CondName;
num_condition = ip.Results.num_condition;
data_all = cell(num_condition,1);
data_tot = cell(num_condition,1);
data_mov_num = cell(num_condition,1);
num_movie = cell(num_condition,1);
Tab = cell(num_condition,1);
for i = 1:num_condition
    if i == 1
        data_all{i} = data1; num_movie{i} = size(data1,2); 
    else 
        data_all{i} = data2; num_movie{i} = size(data2,2);
    end
    data_tot{i} = []; data_mov_num{i} = []; Tab{i} = [];
end

for i = 1:num_condition
    for i_movie = 1:num_movie{i}
    data_tot{i} = cat(1, data_tot{i}, data_all{i}{i_movie});
    num_data = size(data_all{i}{i_movie},1);
    data_mov_num{i} = cat(1, data_mov_num{i}, i_movie*ones(num_data,1));
    end
    Tab{i} = table(data_tot{i},data_mov_num{i});
end


if isempty(BoundedSupport)
    BoundedSupport = false;
end   
if num_condition == 2
fig_pdf = distnObj_plotPDFs_BootCI(Tab{1}, Tab{2}, CondName{1}, CondName{2}, varName, ...
    'BoundedSupport', BoundedSupport, 'save', save);
% StatAnal (2) distnObj_plotCDFs_BootCI 
% pooledTableToDistnObjects
% fprintf('\n\n\n==== Ctrl condition ====\n\n\n')
% DO_1 = pooledTableToDistnObjects(Tab{1});
% fprintf('\n\n\n==== Trt condition ====\n\n\n')
% DO_2 = pooledTableToDistnObjects(Tab{2});
% distnObj_plotCDFs_BootCI
%fig_cdf = distnObj_plotCDFs_BootCI(DO_1, DO_2, CondName{1}, CondName{2}, varName, 'save', save);
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
else
    fig_pdf = distnObj_plotPDFs_BootCI(Tab{1}, [], CondName{1}, [], varName, ...
    'BoundedSupport', BoundedSupport, 'save', save,'trt_exist',false);
% StatAnal (2) distnObj_plotCDFs_BootCI 
% pooledTableToDistnObjects
%fprintf('\n\n\n==== Ctrl condition ====\n\n\n')
%DO_1 = pooledTableToDistnObjects(Tab{1});
% distnObj_plotCDFs_BootCI
%fig_cdf = distnObj_plotCDFs_BootCI(DO_1, [], Name1, [], varName, 'save', save);
end