function [fig_pdf, fig_cdf] = plotPDFconf(data1, data2, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data1', @(x) iscell(x));
ip.addRequired('data2', @(x) iscell(x));
ip.addParameter('BoundedSupport', [4, Inf], @isnumeric); % BoundedSupport = [] (default) or [4, Inf] or [0, Inf]
ip.addParameter('save', false, @islogical);
ip.parse(data1, data2, varargin{:});

save = ip.Results.save;
data1 = ip.Results.data1;
data2 = ip.Results.data2;
BoundedSupport = ip.Results.BoundedSupport;
num_condition = 2;

num_movie = cell(num_condition,1);

num_movie{1} = size(data1,2);
num_movie{2} = size(data2,2);

data1_tot = [];
data2_tot = [];
data1_mov_num = [];
data2_mov_num = [];

for i_movie = 1:num_movie{1}
    data1_tot = cat(1, data1_tot, data1{i_movie});
    num_data = size(data1{i_movie},1);
    data1_mov_num = cat(1, data1_mov_num, i_movie*ones(num_data,1));
end

for i_movie = 1:num_movie{2}
    data2_tot = cat(1, data2_tot, data2{i_movie});
    num_data = size(data2{i_movie},1);
    data2_mov_num = cat(1, data2_mov_num, i_movie*ones(num_data,1));
end


Tab1 = table(data1_tot,data1_mov_num);
Tab2 = table(data2_tot,data2_mov_num);

varName = 'LTmaxIntGrthan150';
Name1 = 'cond1';
Name2 = 'cond2';

if isempty(BoundedSupport)
    BoundedSupport = false;
end   

fig_pdf = distnObj_plotPDFs_BootCI(Tab1, Tab2, Name1, Name2, varName, ...
    'BoundedSupport', BoundedSupport, 'save', save);

% StatAnal (2) distnObj_plotCDFs_BootCI 
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

% pooledTableToDistnObjects
fprintf('\n\n\n==== Ctrl condition ====\n\n\n')
DO_1 = pooledTableToDistnObjects(Tab1);

fprintf('\n\n\n==== Trt condition ====\n\n\n')
DO_2 = pooledTableToDistnObjects(Tab2);

% distnObj_plotCDFs_BootCI


fig_cdf = distnObj_plotCDFs_BootCI(DO_1, DO_2, Name1, Name2, varName, 'save', save);
