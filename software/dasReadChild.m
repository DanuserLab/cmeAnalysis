function [Track_info_child] = dasReadChild(data, master_dir, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('master_dir', @(x) ischar(x));
ip.addParameter('saveTrack_info', true, @islogical);
ip.addParameter('Category', {'Ia'}, @iscell);
ip.parse(data, master_dir, varargin{:});

data = ip.Results.data;
n_ch = size(data(1).channels,2);

master_dir = ip.Results.master_dir;

cd(master_dir);  
% if ip.Results.saveTrack_info == true
%    mkdir('DAS');
% end
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

Category = ip.Results.Category;
%--------------------------------------------------------------------------
movie_num = max(size(data));
%--------------------------------------------------------------------------
fnames = {'Track_num', 'Track_id','LT', 'Max', 'A', 't' , 'x', 'y' ,'frameRate', 'Valid', 'Ib_num', 'sbA', 'ebA', 'sbSigma_r','ebSigma_r','Sigma_r'};

Track_info_child = cell(1,movie_num);
if n_ch > 1
     for i_movie = 1:movie_num
        Track_info_child{i_movie} = cell2struct(cell(size(fnames)), fnames, 2);
     end
end
%==========================================================================
if n_ch > 1
    disp('reading tracks (children) ...');
parfor i_movie = 1:movie_num
    %--------------------------------------------------------------------------
    frame_num = data(i_movie).movieLength;
    track = loadTracks(data(i_movie),'Category', Category);
    track_num = size(track,2);
    nsb = numel(track(1).startBuffer.t);
    neb = numel(track(1).endBuffer.t);
%--------------------------------------------------------------------------
Track_info_child{i_movie}.Track_num = track_num;
Track_info_child{i_movie}.Valid = [];
Track_info_child{i_movie}.LT = zeros(track_num,1);
Track_info_child{i_movie}.Max = [];
Track_info_child{i_movie}.A = nan(track_num,frame_num);
Track_info_child{i_movie}.x = [];
Track_info_child{i_movie}.y = [];
Track_info_child{i_movie}.t = [];
Track_info_child{i_movie}.frameRate = data(i_movie).framerate;
Track_info_child{i_movie}.sbA = nan(track_num,nsb);
Track_info_child{i_movie}.ebA = nan(track_num,neb);
Track_info_child{i_movie}.sbSigma_r = nan(track_num,nsb);
Track_info_child{i_movie}.ebSigma_r = nan(track_num,neb);
Track_info_child{i_movie}.Sigma_r = nan(track_num,frame_num);
%--------------------------------------------------------------------------
    for i_track = 1:track_num
        Track_info_child{i_movie}.LT(i_track) = track(i_track).lifetime_s/Track_info_child{i_movie}.frameRate;
        %Track_info_child{i_movie}.Max(i_track) = max(track(i_track).A(2:n_ch,:));
        Track_info_child{i_movie}.A(i_track,1:Track_info_child{i_movie}.LT(i_track)) = track(i_track).A(2,1:Track_info_child{i_movie}.LT(i_track));
        Track_info_child{i_movie}.Sigma_r(i_track,1:Track_info_child{i_movie}.LT(i_track)) = track(i_track).sigma_r(2,1:Track_info_child{i_movie}.LT(i_track));
        Track_info_child{i_movie}.sbA(i_track,:) = track(i_track).startBuffer.A(2,:);
        Track_info_child{i_movie}.ebA(i_track,:) = track(i_track).endBuffer.A(2,:);
        Track_info_child{i_movie}.sbSigma_r(i_track,:) = track(i_track).startBuffer.sigma_r(2,:);
        Track_info_child{i_movie}.ebSigma_r(i_track,:) = track(i_track).endBuffer.sigma_r(2,:);
    end
    %Track_info_child{i_movie}.A = squeeze(Track_info_child{i_movie}.A);
    %Track_info_child{i_movie}.Sigma_r = squeeze(Track_info_child{i_movie}.Sigma_r);
%--------------------------------------------------------------------------
end
%==========================================================================
end

cd(master_dir); 


