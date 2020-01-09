function [Track_info] = dasRead(data, master_dir, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('master_dir', @(x) ischar(x));
ip.addParameter('saveTrack_info', true, @islogical);
ip.addParameter('Category', {'Ia'}, @iscell);
ip.parse(data, master_dir, varargin{:});

data = ip.Results.data;

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
Track_info = cell(1,movie_num);
fnames = {'Track_num', 'Track_id','LT', 'Max', 'A', 't' , 'x', 'y' ,'frameRate',...
    'Valid', 'Ib_num', 'Id_num', 'sbA', 'ebA', 'sbSigma_r','ebSigma_r','Sigma_r'};

for i_movie = 1:movie_num
    Track_info{i_movie} = cell2struct(cell(size(fnames)), fnames, 2);
end

%==========================================================================
    disp('reading tracks ...');
parfor i_movie = 1:movie_num
    %--------------------------------------------------------------------------
    frame_num = data(i_movie).movieLength;
    track = loadTracks(data(i_movie),'Category', Category);
    track_num = size(track,2);
    track_irreg = loadTracks(data(i_movie),'Category', {'Ib'});
    Track_info{i_movie}.Ib_num = size(track_irreg,2);
    track_irreg = loadTracks(data(i_movie),'Category', {'Id'});
    Track_info{i_movie}.Id_num = size(track_irreg,2);
    nsb = numel(track(1).startBuffer.t);
    neb = numel(track(1).endBuffer.t);
%--------------------------------------------------------------------------
Track_info{i_movie}.Track_num = track_num;
Track_info{i_movie}.Track_id = (1:Track_info{i_movie}.Track_num)';
Track_info{i_movie}.Valid = [];
Track_info{i_movie}.LT = zeros(Track_info{i_movie}.Track_num,1);
Track_info{i_movie}.Max = zeros(Track_info{i_movie}.Track_num,1);
Track_info{i_movie}.A = nan(Track_info{i_movie}.Track_num,frame_num);
Track_info{i_movie}.x = zeros(Track_info{i_movie}.Track_num,frame_num);
Track_info{i_movie}.y = zeros(Track_info{i_movie}.Track_num,frame_num);
Track_info{i_movie}.t = zeros(Track_info{i_movie}.Track_num,frame_num);
Track_info{i_movie}.frameRate = data(i_movie).framerate;
Track_info{i_movie}.sbA = nan(Track_info{i_movie}.Track_num,nsb);
Track_info{i_movie}.ebA = nan(Track_info{i_movie}.Track_num,neb);
Track_info{i_movie}.sbSigma_r = nan(Track_info{i_movie}.Track_num,nsb);
Track_info{i_movie}.ebSigma_r = nan(Track_info{i_movie}.Track_num,neb);
Track_info{i_movie}.Sigma_r = nan(Track_info{i_movie}.Track_num,frame_num);
%--------------------------------------------------------------------------
i_ch = 1;
    for i_track = 1:Track_info{i_movie}.Track_num
        Track_info{i_movie}.LT(i_track) = track(i_track).lifetime_s/Track_info{i_movie}.frameRate;
        Track_info{i_movie}.Max(i_track) = max(track(i_track).A(i_ch,:));
        Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)) = track(i_track).A(i_ch,1:Track_info{i_movie}.LT(i_track));
        Track_info{i_movie}.Sigma_r(i_track,1:Track_info{i_movie}.LT(i_track)) = track(i_track).sigma_r(i_ch,1:Track_info{i_movie}.LT(i_track));
        %Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)) = track(i_track).A(1:Track_info{i_movie}.LT(i_track))./exp(-0.08*track(i_track).t(1:Track_info{i_movie}.LT(i_track)));
        %Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)) = track(i_track).A(Track_info{i_movie}.LT(i_track):-1:1);
        Track_info{i_movie}.x(i_track,1:Track_info{i_movie}.LT(i_track)) = floor(track(i_track).x(i_ch,:)+0.5);
        Track_info{i_movie}.y(i_track,1:Track_info{i_movie}.LT(i_track)) = floor(track(i_track).y(i_ch,:)+0.5);
        Track_info{i_movie}.t(i_track,1:Track_info{i_movie}.LT(i_track)) = track(i_track).t(:);
        %Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)) = Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)) ./ exp(-0.0068*Track_info{i_movie}.t(i_track,1:Track_info{i_movie}.LT(i_track)))/611*109;
        %Track_info{i_movie}.Max(i_track) = max(Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track)));
        Track_info{i_movie}.sbA(i_track,:) = track(i_track).startBuffer.A(i_ch,:);
        Track_info{i_movie}.ebA(i_track,:) = track(i_track).endBuffer.A(i_ch,:);
        Track_info{i_movie}.sbSigma_r(i_track,:) = track(i_track).startBuffer.sigma_r(i_ch,:);
        Track_info{i_movie}.ebSigma_r(i_track,:) = track(i_track).endBuffer.sigma_r(i_ch,:);
    end
    %--------------------------------------------------------------------------
%==========================================================================

%==========================================================================
end




cd(master_dir); 


