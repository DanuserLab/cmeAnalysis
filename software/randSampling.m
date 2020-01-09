

function  data_samp = randSampling(data, ratio, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isnumeric(x));
ip.addRequired('ratio', @(x) isnumeric(x));
ip.addParameter('dim', 1, @isnumeric);
ip.parse(data, ratio, varargin{:});


data = ip.Results.data;

ratio = ip.Results.ratio;

dim = ip.Results.dim;

n = size(data);

n_to_samp = n(dim);

samp_num = floor(ratio*n_to_samp+0.5);

       [~,id_temp] = sort(rand(n_to_samp,1));
       id_temp(samp_num+1:n_to_samp) = [];
       data_samp = data(id_temp,:);





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
