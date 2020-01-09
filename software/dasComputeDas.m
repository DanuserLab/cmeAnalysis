
function DAS_all = dasComputeDas(data, Track_info, DAS_all, Dfunc, factor, A_max_all, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addRequired('DAS_all', @(x) isstruct(x));
ip.addRequired('Dfunc', @(x) ismatrix(x));
ip.addRequired('factor', @(x) isnumeric(x) );
ip.addRequired('A_max_all', @(x) isnumeric(x) );
ip.addParameter('bin_A', 1, @isnumeric);
ip.addParameter('calculate_area', false, @islogical);

ip.addParameter('MaskPath', ['Detection' filesep 'cellmask.tif'], @ischar);
ip.parse(data, Track_info, DAS_all, Dfunc, factor, A_max_all,varargin{:});

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
A_max_all = ip.Results.A_max_all;
bin_A = ip.Results.bin_A;
A_bin_num = A_max_all/bin_A;
MaskPath = ['Detection' filesep 'cellmask.tif'];
DAS_all = ip.Results.DAS_all;
factor = ip.Results.factor;
calculate_area = ip.Results.calculate_area;
%--------------------------------------------------------------------------------------------------------
disp('computing DAS ...');
Track_info = ip.Results.Track_info;
Dfunc_temp = cell(movie_num,1);

for i_movie = 1:movie_num
    Dfunc_temp{i_movie} = ip.Results.Dfunc;
end

if calculate_area == true
cellAreaTime = zeros(movie_num, 0);
for i_movie = 1:movie_num
    %---------------------------------
    px = data(i_movie).pixelSize / data(i_movie).M; % pixels size in object space
    mask = logical(readtiff([data(i_movie).source MaskPath]));
    cellAreaTime(i_movie) = sum(mask(:)) * px^2 * 1e12 * data(i_movie).movieLength/60; % in �m^2
    %---------------------------------
end
end

DAS = cell(movie_num,1);
DAS_var = cell(movie_num,1);
DAS_2 = cell(movie_num,1);
DAS_3 = cell(movie_num,1);
LT_thre = 1;
frame_start = 1;
for i_movie = 1:movie_num

   DAS_temp = zeros(Track_info{i_movie}.Track_num,data(i_movie).movieLength);
   for i_track = 1:Track_info{i_movie}.Track_num
       A_int = floor(Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track))*factor(1)+factor(2) / bin_A +0.5);    

       %frame_start = Track_info{i_movie}.LT(i_track)-LT_thre;
       for i_frame = frame_start: Track_info{i_movie}.LT(i_track)
           if (A_int(i_frame) > 0) && (A_int(i_frame) < A_bin_num)   
            DAS_temp(i_track,i_frame) = Dfunc_temp{i_movie}(A_int(i_frame),i_frame);  
           end
       end
       

   end
   
   %DAS{i_movie} = DAS_temp/LT_thre;
    DAS{i_movie} = sum(DAS_temp,2)./(Track_info{i_movie}.LT-1);
    %DAS_var{i_movie} = sqrt((DAS_temp - DAS{i_movie}).^2./DAS{i_movie}.^2);
    %DAS_var{i_movie} = sum(DAS_var{i_movie},2)./(Track_info{i_movie}.LT-1);
    DAS_var{i_movie} = zeros(Track_info{i_movie}.Track_num,1);
    DAS_2{i_movie} = zeros(Track_info{i_movie}.Track_num,1);
    DAS_3{i_movie} = zeros(Track_info{i_movie}.Track_num,1);
    for i_track = 1:Track_info{i_movie}.Track_num
    DAS_var{i_movie}(i_track) = log((max(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)))-min(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track))))/Track_info{i_movie}.LT(i_track));
    %DAS_var{i_movie}(i_track) = var(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)));
    %[M,I] = max(diff(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track))));
    %DAS_test{i_movie}(i_track) = I;
    %DAS_test{i_movie}(i_track) = min(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)));
    DAS_2{i_movie}(i_track) = mean(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)).^2);
    DAS_3{i_movie}(i_track) = mean(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)).^3);
    %DAS_var{i_movie}(i_track) = quantile(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)),0.90) - quantile(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)),0.10);
    %DAS_var{i_movie}(i_track) = var(DAS_temp(i_track,1:Track_info{i_movie}.LT(i_track)));
    end
end
DAS_all.DAS = [];
DAS_all.DAS_var = [];
DAS_all.DAS_2 = [];
DAS_all.DAS_3 = [];
for i_movie = 1:movie_num
    DAS_all.DAS = cat(1,DAS_all.DAS,DAS{i_movie});
    DAS_all.DAS_var = cat(1,DAS_all.DAS_var,DAS_var{i_movie});
    DAS_all.DAS_2 = cat(1,DAS_all.DAS_2,DAS_2{i_movie});
    DAS_all.DAS_3 = cat(1,DAS_all.DAS_3,DAS_3{i_movie});
end
DAS_all.ImaxAll = A_max_all;
DAS_all.DAS_var(isinf(DAS_all.DAS_var)) = log(0.00000000000000000000000000000000000000000001);



