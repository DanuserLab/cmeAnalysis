function [data_all] = poolExperiments(varargin)
% this function pools different experiment structures into one, e.g.
% different control datasets of different origins
% INPUT:    various input experiments, e.g. data1, data2, data3
%           where each experiment has multiple fields (e.g. .source,
%           .framerate, .movieLength)
% OUTPUT:   merged experiment file
%
% Dinah Loerke, March 15, 2009
% Last Modified: Francois Aguet, 01/21/2010
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

for i=1:length(varargin)
    
    data_curr = varargin{i};
    
    if i==1
        names = fieldnames(data_curr);
        data_all = data_curr;
    else        
        clen = length(data_all);
        for n = 1:length(data_curr)
            curr_struct = data_curr(n);
            for k = 1:length(names)
                data_all    = setfield(data_all, {1,clen+n}, names{k}, curr_struct.(names{k}));
            end
        end
    end
end