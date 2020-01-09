

function rand_track = plotExampleTrack(data, Track_info,id, fig, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addRequired('id', @(x) isnumeric(x));
ip.addRequired('fig', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('plot_all', false, @islogical);
ip.addParameter('plot_frame', false, @islogical);
ip.addParameter('i_movie', 1, @isnumeric);
ip.addParameter('montage_x_num', 15, @isnumeric);
ip.addParameter('resize_factor', 1, @isnumeric);
ip.addParameter('plot_partial', false, @islogical);
ip.addParameter('subplot_id', [1,2,1:2], @isnumeric);
ip.addOptional('track_id', [], @isnumeric);
ip.parse(data, Track_info,id, fig,pm,varargin{:});

data = ip.Results.data;
Track_info = ip.Results.Track_info;
resize_factor = ip.Results.resize_factor;
montage_x_num = ip.Results.montage_x_num;
id = ip.Results.id;
fig = ip.Results.fig;
pm = ip.Results.pm;
subplot_id = ip.Results.subplot_id;
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
i_movie = ip.Results.i_movie;
bandsize = 0;
%--------------------------------------------------------------------------------------------------------
image_name = data(i_movie).framePaths;
frame_tot = data(i_movie).movieLength;

%--------------------------------------------------------------------------------------------------------
patch_size = 15;
patch_size_half = (patch_size-1)/2;
%--------------------------------------------------------------------------------------------------------
clims = [150 300];
%--------------------------------------------------------------------------------------------------------
exp_num = max(size(id));

rand_track = struct('id', [],'LT', [], 't', [], 'A', [], 'x', [], 'y', [],'image_map', []);

rand_track.id = Track_info{i_movie}.Track_id(id);
rand_track.t = cell(exp_num,1);
rand_track.x = cell(exp_num,1);
rand_track.y = cell(exp_num,1);
rand_track.image_map = cell(exp_num,1);

%--------------------------------------------------------------------------------------------------------
frame_all = [];
for i_exp = 1:exp_num
    i_track = rand_track.id(i_exp);
    rand_track.LT(i_exp) = Track_info{i_movie}.LT(i_track);
    rand_track.t{i_exp} = Track_info{i_movie}.t(i_track,1:Track_info{i_movie}.LT(i_track));
    rand_track.A{i_exp} = Track_info{i_movie}.A(i_track,1:Track_info{i_movie}.LT(i_track));
    rand_track.x{i_exp} = Track_info{i_movie}.x(i_track,1:Track_info{i_movie}.LT(i_track));
    rand_track.y{i_exp} = Track_info{i_movie}.y(i_track,1:Track_info{i_movie}.LT(i_track));
    rand_track.image_map{i_exp} = cell(rand_track.LT(i_exp),1);
    frame_all = cat(1, frame_all, rand_track.t{i_exp}');
end

if pm.fig1_mosaic
for i_frame = 1:frame_tot
   if ~isempty(frame_all(frame_all == i_frame))
        image_map = imread(image_name{1},i_frame);
        image_map = image_map';
        image_map_mirror = image_map;

        for i_exp = 1:exp_num           
            for i_t = 1:rand_track.LT(i_exp)
                if rand_track.t{i_exp}(i_t) == i_frame
                x_temp = rand_track.x{i_exp}(i_t)*resize_factor;
                y_temp = rand_track.y{i_exp}(i_t)*resize_factor;
                rand_track.image_map{i_exp}{i_t} = image_map(x_temp-patch_size_half:x_temp+patch_size_half,y_temp-patch_size_half:y_temp+patch_size_half);
                for i_band = -bandsize:bandsize
                image_map_mirror(x_temp+patch_size_half+i_band,y_temp-patch_size_half:y_temp+patch_size_half) = clims(2);
                image_map_mirror(x_temp-patch_size_half+i_band,y_temp-patch_size_half:y_temp+patch_size_half) = clims(2);
                image_map_mirror(x_temp-patch_size_half:x_temp+patch_size_half,y_temp-patch_size_half+i_band) = clims(2);
                image_map_mirror(x_temp-patch_size_half:x_temp+patch_size_half,y_temp+patch_size_half+i_band) = clims(2);
                end

                %image_map_mirror(x_temp,y_temp) = clims(2);
                end
            end
        end
     if ip.Results.plot_frame == true
     colormap('hot')
     imagesc(image_map_mirror,clims);
     pause();
     end
   end
end
end

if pm.fig1_mosaic
    montage = [];
    for i_exp = 1:exp_num  
    montage_x = [];
        for i_t = 1:montage_x_num
            if i_t <= rand_track.LT(i_exp)
                montage_x = cat(2,montage_x,rand_track.image_map{i_exp}{i_t});
            else
                montage_x = cat(2,montage_x,uint16(zeros(patch_size,patch_size)));
            end
        end
     montage = cat(1, montage, montage_x);
    end
     figure(fig{1});

     subplot(subplot_id(1),subplot_id(2),subplot_id(3))
     colormap('hot')
     imagesc(montage,clims);
     axis off
     
     subplot(subplot_id(1),subplot_id(2),subplot_id(4))
     hold on
     for i_exp = 1:exp_num  
         plot(rand_track.A{i_exp}(:),'Linewidth', 0.5)
     end
     xlabel('$t$ (s)','interpreter','latex')
     ylabel('$I$ (a.u.)','interpreter','latex')
     ylim([0 300])
     xlim([0 200])
     hold off
else

     figure(fig{1});
     hold on
     for i_exp = 1:exp_num  
         plot(rand_track.A{i_exp}(:),'Linewidth', 0.5)
     end
     xlabel('$t$ (s)','interpreter','latex')
     ylabel('$I$ (a.u.)','interpreter','latex')
     ylim([0 300])
     xlim([0 200])

     hold off
end

%      length_x = 0.2+0.8/20*montage_x_num;
%      length_y = 1.0;
%      set(gcf, 'Units', 'Normalized', 'OuterPosition', [0., 0., length_x, length_y]);
%      set(gca,'fontsize',17, 'FontWeight', 'bold')
%      set(gca,'Linewidth',2)





