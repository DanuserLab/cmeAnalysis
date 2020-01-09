
function DAS_all = dasDataProcess(data, Track_info, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addParameter('MaskPath', ['Detection' filesep 'cellmask.tif'], @ischar);
ip.addParameter('A10_norm', 60, @isnumeric); %CLC-eGFP: 60, alpha-eGFP: 40
ip.addParameter('LT_thre', 5, @isnumeric);
ip.parse(data, Track_info, varargin{:});

%--------------------------------------------------------------------------
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
movie_num = size(ip.Results.data,2);
data = ip.Results.data;
Track_info_temp = ip.Results.Track_info;
MaskPath = ['Detection' filesep 'cellmask.tif'];
A10_norm = ip.Results.A10_norm;
LT_thre = ip.Results.LT_thre;
%--------------------------------------------------------------------------------------------------------
cellAreaTime = zeros(movie_num, 0);
for i_movie = 1:movie_num
    %---------------------------------
    px = data(i_movie).pixelSize / data(i_movie).M; % pixels size in object space
    mask = logical(readtiff([data(i_movie).source MaskPath]));
    cellAreaTime(i_movie) = sum(mask(:)) * px^2 * 1e12 * data(i_movie).movieLength/60; % in �m^2
    %---------------------------------
end

DAS_all = struct('DAS',[], 'LT', [], 'MaxI', [],'MovieNum', [], 'TrackID', [],...
    'cellAreaTime', [], 'DAS_var', [], 'DAS_2', [],'DAS_3', [],'A5', [], ...
    'factor', [], 'Ib_num', [], 'Id_num', [], 'ImaxAll',[]);

for i_movie = 1:movie_num
    DAS_all.LT = cat(1,DAS_all.LT,Track_info_temp{i_movie}.LT(:));
    DAS_all.MaxI = cat(1,DAS_all.MaxI,Track_info_temp{i_movie}.Max(:));
    DAS_all.MovieNum = cat(1,DAS_all.MovieNum,ones(Track_info_temp{i_movie}.Track_num,1)*i_movie);
    DAS_all.TrackID = cat(1,DAS_all.TrackID,Track_info_temp{i_movie}.Track_id(:));
    DAS_all.cellAreaTime = cat(1,DAS_all.cellAreaTime,cellAreaTime(i_movie));
    DAS_all.A5 = cat(1, DAS_all.A5, [Track_info_temp{i_movie}.A(:,4) Track_info_temp{i_movie}.A(:,5)]);
    DAS_all.Ib_num = cat(1, DAS_all.Ib_num, Track_info_temp{i_movie}.Ib_num);
    DAS_all.Id_num = cat(1, DAS_all.Id_num, Track_info_temp{i_movie}.Id_num);
end
DAS_all.DAS = [];



