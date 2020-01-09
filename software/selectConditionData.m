%[data] = selectConditionData(varargin) displays a GUI for selection of individual movies from a list
% If no input is provided, this function first loads data sets using loadConditionData.
%
% Input (optional)
%         data : movie data structure returned by loadConditionData()
%
% Output
%         data : selected subset of movie data structures
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

% Daniel Nunez, 2008 (last modified 06/05/2013 by F.A.)

function [data] = selectConditionData(data)

% Load data sets
if nargin<1 || isempty(data)
    data = loadConditionData();
    loadData = true;
    while loadData
        str = [];
        while ~any(strcmpi(str, {'y','n'}))
            str = input('Load additional ''condition'' folders? [y/n]: ', 's');
        end
        if strcmpi(str, 'y')
            [data] = [data loadConditionData()];
        else
            loadData = false;
        end
    end
end

idx = listSelectGUI({data.source}, [], 'move', 1);
data = data(idx);
