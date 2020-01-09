

function [fig_exist] = dasCohorts(data, DAS_all, idx,dir_alt, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @iscell);
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('idx', @(x) iscell(x));
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('ShowVariation', true, @islogical);
ip.addParameter('FillMode', 'SEM', @(x) any(strcmpi(x, {'SEM', 'pct'})));
ip.addParameter('FrontLayer', false, @islogical);
ip.addParameter('ShowBackground', false, @islogical);
ip.addParameter('Rescale', true, @islogical);
ip.addParameter('RescalingReference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
ip.addParameter('ScaleSlaveChannel', false, @islogical);
ip.addParameter('MaxIntensityThreshold', 0);
ip.addParameter('ExcludeVisitors', true, @islogical);
ip.addParameter('SlaveName', [], @iscell);
ip.addParameter('ChannelNames', []);
ip.addParameter('LineStyle', '-');
ip.addParameter('Hues', []);
ip.addParameter('Colormap', []);
ip.addParameter('ColormapFill', []);
ip.addParameter('DisplayMode', 'screen');
ip.addParameter('DisplayAll', false);
ip.addParameter('TrackIndex', []);
ip.addParameter('Cutoff_f', 5);
ip.addParameter('Alpha', 0.05);
ip.addParameter('YTick', []);
ip.addParameter('YLim', []);
ip.addParameter('Parent', []);
ip.addParameter('RemoveOutliers', false, @islogical);
ip.addParameter('PlotXLabel', true, @islogical);
ip.addParameter('ShowPct', true, @islogical);
ip.addParameter('ShowStats', false, @islogical);
ip.addParameter('AvgFun', @nanmean, @(x) isa(x, 'function_handle'));
ip.addParameter('LftDataName', 'lifetimeData.mat');
ip.addParameter('AmplitudeCorrection', []);
ip.addParameter('Align', 'left', @(x) any(strcmpi(x, {'right', 'left'})));
ip.addParameter('fig_exist', [], @iscell);
ip.addParameter('re_read', true,@islogical);
ip.addParameter('b', 5, @isnumeric);
ip.addParameter('master_adj', 0, @isnumeric);
%==========================================================================
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
ip.parse(data, DAS_all, idx, dir_alt, pm, varargin{:});
%==========================================================================
data = ip.Results.data;
DAS_all = ip.Results.DAS_all;
idx = ip.Results.idx;
dir_alt = ip.Results.dir_alt;
fig_exist = ip.Results.fig_exist;
plot_visitor = pm.plot_visitor;
pm = ip.Results.pm;
MaxIntensityThreshold = ip.Results.MaxIntensityThreshold;
master_adj = ip.Results.master_adj;
num_clus = 3;
ScaleSlaveChannel = ip.Results.ScaleSlaveChannel;
if pm.cohort_diff == false
    re_read = ip.Results.re_read;
else
    re_read = true;
end

cohortBounds = cell(num_clus,1);
cohortBounds{1} = pm.CohortBounds_s_ccp;
cohortBounds{2} = pm.CohortBounds_s_visitor;
cohortBounds{3} = pm.CohortBounds_s_ftn;

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

b = ip.Results.b;
framerate = data{1}(1).framerate;

% # data points in cohort (including buffer frames)
iLength = cell(num_clus);
cT = cell(num_clus);
nc = cell(num_clus);
for i_c = 1:num_clus
    nc{i_c} = numel(cohortBounds{i_c})-1;
iLength{i_c} = arrayfun(@(c) floor(mean(cohortBounds{i_c}([c c+1]))/framerate) + 2*b, 1:nc{i_c});
% time vectors for cohorts
if strcmpi(ip.Results.Align, 'left')
    cT{i_c} = arrayfun(@(i) (-b:i-b-1)*framerate, iLength{i_c}, 'unif', 0);
else
    cT{i_c} = arrayfun(@(i) (-i+b+1:b)*framerate, iLength{i_c}, 'unif', 0);
end
end


nCh = max(size(data{1}(1).channels));
num_condition = max(size(data));


%==========================================================================

%==========================================================================
if (~(exist(fullfile(dir_alt, 'DAS_cohort.mat'), 'file') == 2)) || (re_read == true)
%==========================================================================
disp('reading tracks...');
S = cell(num_condition,1);
for i_condition = 1:num_condition
    if i_condition == 1
        file_name_temp = [dir_alt filesep 'ctrl' filesep 'Track_info.mat'];
    elseif i_condition==2
        file_name_temp = [dir_alt filesep 'treated' filesep 'Track_info.mat'];
    else
        file_name_temp = [dir_alt filesep 'treated' num2str(i_condition) filesep 'Track_info.mat'];
    end
    S{i_condition} = load(file_name_temp);
end
disp('finished reading tracks');
clear i_condition;
clear file_name_temp;

if nCh > 1
disp('reading tracks...');
S_c = cell(num_condition,1);
for i_condition = 1:num_condition
    if i_condition == 1
        file_name_temp = [dir_alt filesep 'ctrl' filesep 'Track_info_child.mat'];
    elseif i_condition==2
        file_name_temp = [dir_alt filesep 'treated' filesep 'Track_info_child.mat'];
    else
        file_name_temp = [dir_alt filesep 'treated' num2str(i_condition) filesep 'Track_info_child.mat'];
    end
    S_c{i_condition} = load(file_name_temp);
end
disp('finished reading tracks');
clear i_condition;
clear file_name_temp;
end
%==========================================================================
res = cell(num_condition,num_clus);

for i_condition = 1:num_condition
    for i_c = 1:num_clus
    res{i_condition,i_c} = struct(...
    'A', [], ...
    'x', [], ...
    'y', [], ...
    'sbA', [],...
    'ebA', [],...
    'sbSigma_r', [],...
    'ebSigma_r', [],...
    'Sigma_r', [],...
    'interpTracks', [],...
    'interpSigLevel', [],...
    'n_movie', [],...
    'A_mean_t', [],...
    'A_mean_err', [],...
    'LT',[]...
     );
    end
end
if nCh > 1
    res_c = res;
else
    res_c = [];
end
%==========================================================================
for i_condition = 1:num_condition
    %-----------------------------------------------------------------------
    Track_info = S{i_condition}.Track_info;
    %-----------------------------------------------------------------------
    for i_c = 1:num_clus
    
    res{i_condition,i_c}.A = [];
    res{i_condition,i_c}.LT = [];

    n_movie = max(size(Track_info));
    res{i_condition,i_c}.n_movie = n_movie;
    %----------------------------------------------------------------------
        res{i_condition,i_c}.A = cell(n_movie,1);
        res{i_condition,i_c}.Sigma_r = cell(n_movie,1);
        res{i_condition,i_c}.LT = cell(n_movie,1);
        res{i_condition,i_c}.x = cell(n_movie,1);
        res{i_condition,i_c}.y = cell(n_movie,1);
        res{i_condition,i_c}.sbA = cell(n_movie,1);
        res{i_condition,i_c}.ebA = cell(n_movie,1);
        res{i_condition,i_c}.sbSigma_r = cell(n_movie,1);
        res{i_condition,i_c}.ebSigma_r = cell(n_movie,1);
    %----------------------------------------------------------------------
    for i_mov = 1:n_movie
        id_temp = (DAS_all{i_condition}.MovieNum == i_mov) & (idx{i_condition} == i_c)...
                  & (DAS_all{i_condition}.MaxI > MaxIntensityThreshold);
        res{i_condition,i_c}.A{i_mov} = Track_info{i_mov}.A(DAS_all{i_condition}.TrackID(id_temp),:);
        res{i_condition,i_c}.Sigma_r{i_mov} = Track_info{i_mov}.Sigma_r(DAS_all{i_condition}.TrackID(id_temp),:);
        res{i_condition,i_c}.LT{i_mov} = DAS_all{i_condition}.LT(id_temp);
        res{i_condition,i_c}.x{i_mov} = mean(Track_info{i_mov}.x(DAS_all{i_condition}.TrackID(id_temp),:),2,'omitnan');
        res{i_condition,i_c}.y{i_mov} = mean(Track_info{i_mov}.y(DAS_all{i_condition}.TrackID(id_temp),:),2,'omitnan');
        res{i_condition,i_c}.sbA{i_mov} = Track_info{i_mov}.sbA(DAS_all{i_condition}.TrackID(id_temp),:);
        res{i_condition,i_c}.ebA{i_mov} = Track_info{i_mov}.ebA(DAS_all{i_condition}.TrackID(id_temp),:);
        res{i_condition,i_c}.sbSigma_r{i_mov} = Track_info{i_mov}.sbSigma_r(DAS_all{i_condition}.TrackID(id_temp),:);
        res{i_condition,i_c}.ebSigma_r{i_mov} = Track_info{i_mov}.ebSigma_r(DAS_all{i_condition}.TrackID(id_temp),:);           
    end
    end
    %-----------------------------------------------------------------------	
end
if nCh > 1
for i_condition = 1:num_condition
    %-----------------------------------------------------------------------
    Track_info = S_c{i_condition}.Track_info_child;
    %-----------------------------------------------------------------------
    for i_c = 1:num_clus
    
    res_c{i_condition,i_c}.A = [];
    res_c{i_condition,i_c}.LT = [];

    n_movie = max(size(Track_info));
    res_c{i_condition,i_c}.n_movie = n_movie;
    %----------------------------------------------------------------------
        res_c{i_condition,i_c}.A = cell(n_movie,1);
        res_c{i_condition,i_c}.Sigma_r = cell(n_movie,1);
        res_c{i_condition,i_c}.LT = cell(n_movie,1);
        res_c{i_condition,i_c}.x = cell(n_movie,1);
        res_c{i_condition,i_c}.y = cell(n_movie,1);
        res_c{i_condition,i_c}.sbA = cell(n_movie,1);
        res_c{i_condition,i_c}.ebA = cell(n_movie,1);
        res_c{i_condition,i_c}.sbSigma_r = cell(n_movie,1);
        res_c{i_condition,i_c}.ebSigma_r = cell(n_movie,1);
    %----------------------------------------------------------------------
    for i_mov = 1:n_movie

        id_temp = (DAS_all{i_condition}.MovieNum == i_mov) & (idx{i_condition} == i_c) ...
                  & (DAS_all{i_condition}.MaxI > MaxIntensityThreshold);
        res_c{i_condition,i_c}.A{i_mov} = Track_info{i_mov}.A(DAS_all{i_condition}.TrackID(id_temp),:);
        res_c{i_condition,i_c}.Sigma_r{i_mov} = Track_info{i_mov}.Sigma_r(DAS_all{i_condition}.TrackID(id_temp),:);
        res_c{i_condition,i_c}.LT{i_mov} = DAS_all{i_condition}.LT(id_temp);
        res_c{i_condition,i_c}.sbA{i_mov} = Track_info{i_mov}.sbA(DAS_all{i_condition}.TrackID(id_temp),:);
        res_c{i_condition,i_c}.ebA{i_mov} = Track_info{i_mov}.ebA(DAS_all{i_condition}.TrackID(id_temp),:);
        res_c{i_condition,i_c}.sbSigma_r{i_mov} = Track_info{i_mov}.sbSigma_r(DAS_all{i_condition}.TrackID(id_temp),:);
        res_c{i_condition,i_c}.ebSigma_r{i_mov} = Track_info{i_mov}.ebSigma_r(DAS_all{i_condition}.TrackID(id_temp),:);           
    end
    end
    %-----------------------------------------------------------------------
end
end
clear S;
clear S_c;
%==========================================================================
%==========================================================================
%==========================================================================
% loop through data sets, generate cohorts for each
for i = 1:num_condition
    for i_c = 1:num_clus
        %------------------------------------------------------------------
        res{i,i_c}.interpTracks = cell(res{i,i_c}.n_movie,nc{i_c});
        res{i,i_c}.interpSigLevel = cell(res{i,i_c}.n_movie,nc{i_c});
        for i_mov = 1: res{i,i_c}.n_movie
        % interpolate tracks to mean cohort length
        for c = 1:nc{i_c} % cohorts
            cidx = find(cohortBounds{i_c}(c)<=res{i,i_c}.LT{i_mov}/framerate & res{i,i_c}.LT{i_mov}/framerate<cohortBounds{i_c}(c+1) );           
            nt = numel(cidx);
            if nt>0
                interpTracks = zeros(nt,iLength{i_c}(c));
                sigma_rMat = zeros(nt,iLength{i_c}(c));
                cLengths = res{i,i_c}.LT{i_mov}(cidx);
                % loop through track lengths within cohort
                for t = 1:nt                  
                    A = [res{i,i_c}.sbA{i_mov}(cidx(t),:) res{i,i_c}.A{i_mov}(cidx(t),1:cLengths(t)) res{i,i_c}.ebA{i_mov}(cidx(t),:)];
                    bgr = [res{i,i_c}.sbSigma_r{i_mov}(cidx(t),:) res{i,i_c}.Sigma_r{i_mov}(cidx(t),1:cLengths(t)) res{i,i_c}.ebSigma_r{i_mov}(cidx(t),:)];
                    if pm.cohort_diff == true
                     A = diff(A,[],2);                     
                    end
                    A(isnan(A)) = nanmin(A);
                    bgr(isnan(bgr)) = nanmin(bgr);                   
                    % interpolate to mean length
                    xi = linspace(1,cLengths(t)+2*b, iLength{i_c}(c));
                    interpTracks(t,:) = binterp(A, xi);
                    %interpTracks(t,:) = normalize(interpTracks(t,:));
                    sigma_rMat(t,:) = binterp(bgr, xi);
                end                
                res{i,i_c}.interpTracks{i_mov,c} = interpTracks;
                res{i,i_c}.interpSigLevel{i_mov,c} = kLevel*sigma_rMat;               
            else
                res{i,i_c}.interpTracks{i_mov,c} = NaN(1,iLength{i_c}(c));
                res{i,i_c}.interpSigLevel{i_mov,c} = NaN(1,iLength{i_c}(c));
            end
        end
        end
        %------------------------------------------------------------------
        if nCh > 1
        %------------------------------------------------------------------
        res_c{i,i_c}.interpTracks = cell(res{i,i_c}.n_movie,nc{i_c});
        res_c{i,i_c}.interpSigLevel = cell(res_c{i,i_c}.n_movie,nc{i_c});
        for i_mov = 1: res_c{i,i_c}.n_movie
        % interpolate tracks to mean cohort length
        for c = 1:nc{i_c} % cohorts
            cidx = find(cohortBounds{i_c}(c)<=res_c{i,i_c}.LT{i_mov}/framerate & res_c{i,i_c}.LT{i_mov}/framerate<cohortBounds{i_c}(c+1) );           
            nt = numel(cidx);
            if nt>0
                interpTracks = zeros(nt,iLength{i_c}(c));
                sigma_rMat = zeros(nt,iLength{i_c}(c));
                cLengths = res_c{i,i_c}.LT{i_mov}(cidx);
                % loop through track lengths within cohort
                for t = 1:nt
                    A = [res_c{i,i_c}.sbA{i_mov}(cidx(t),:) res_c{i,i_c}.A{i_mov}(cidx(t),1:cLengths(t)) res_c{i,i_c}.ebA{i_mov}(cidx(t),:)];
                    bgr = [res_c{i,i_c}.sbSigma_r{i_mov}(cidx(t),:) res_c{i,i_c}.Sigma_r{i_mov}(cidx(t),1:cLengths(t)) res_c{i,i_c}.ebSigma_r{i_mov}(cidx(t),:)];
                    if pm.cohort_diff == true
                     A = diff(A,[],2);                     
                    end
                    A(isnan(A)) = nanmin(A);
                    bgr(isnan(bgr)) = nanmin(bgr);                   
                    % interpolate to mean length
                    xi = linspace(1,cLengths(t)+2*b, iLength{i_c}(c));
                    interpTracks(t,:) = binterp(A, xi);
                    %interpTracks(t,:) = normalize(interpTracks(t,:));
                    sigma_rMat(t,:) = binterp(bgr, xi);
                end                
                res_c{i,i_c}.interpTracks{i_mov,c} = interpTracks;
                res_c{i,i_c}.interpSigLevel{i_mov,c} = kLevel*sigma_rMat;               
            else
                res_c{i,i_c}.interpTracks{i_mov,c} = NaN(1,iLength{i_c}(c));
                res_c{i,i_c}.interpSigLevel{i_mov,c} = NaN(1,iLength{i_c}(c));
            end
        end
        
        end
        %------------------------------------------------------------------
        end
    end
end
%==========================================================================
%==========================================================================
%==========================================================================
if pm.cohort_diff == false
save([dir_alt filesep 'DAS_cohort.mat'],'res', 'res_c','-v7.3');
end
%==========================================================================
else
    S = load([dir_alt filesep 'DAS_cohort.mat']);
    res = S.res;
    res_c = S.res_c;
    clear S;
end
%==========================================================================


% Set colormap depending on # channels
cmap = pm.color_clus;
cv = cmap;


%==================================================
% Plot cohorts
%==================================================

A = cell(nCh,nc{i_c});
if isempty(fig_exist)
        fig_exist = cell(num_condition, 1);
        for i = 1:num_condition           
            fig_exist{i} = figure; 
        end
end
if pm.cohort_norm
    ha = cell(num_condition, 2);
else
    ha = cell(num_condition, 1);
end
        for i = 1:num_condition

                fig_tem=figure(fig_exist{i});hold on;
set(fig_tem,'PaperUnits','centimeters');
set(fig_tem,'PaperPosition',[0 0 5 5]);
set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');

            if ~pm.cohort_norm
                ha{i} = gca;
                xlabel('t (s)'); ylabel([pm.ch1_name ' Int. (a.u.)']); 
                set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
                if ~pm.coh_buff 
                xlim([0 cT{1}{end}(end-b+1)]);
                end
                ylim([10 70]);
            else
                ha{i,1} = gca;
                %xlabel('t (s)'); ylabel([pm.ch1_name ' Int. (a.u.)']); 
                set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
                if ~pm.coh_buff 
                xlim([0 cT{1}{end}(end-b+1)]);
                end
                ylim([10 70]);
                ax1 = gca;
                ha{i,2} = axes('Position',ax1.Position,'Linewidth',1,...
                            'XAxisLocation','bottom',...
                            'YAxisLocation','right',...
                            'Color','none');
                %ylabel([pm.ch2_name ' Int. (a.u.)']); 
                set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
                %set(gca,'xticklabel',[]);
                if ~pm.coh_buff 
                xlim([0 cT{1}{end}(end-b+1)]);
                end
                ylim([5 35]);
            end
            
            figure(fig_exist{i}); 

        end




sf = ones(1,nCh);
for i = 1:num_condition
    d_i_c = 1;
    if (plot_visitor == false)
        d_i_c = 2;
    end
    for i_c = 1:d_i_c:num_clus 
        %-----------------------
        if strcmpi(ip.Results.Align, 'left')
             tOffset = @(c) 0;
        else
             tOffset = @(c) cT{i_c}{c}(end)-b*framerate;
        end
        %-----------------------
        i_ch = 1;
        
           res{i,i_c}.A_mean_t = cell(nc{i_c},1); 
           res{i,i_c}.A_mean_err = cell(nc{i_c},1);
        for c = 1:nc{i_c}
                AMat = [];
                for i_mov = 1: res{i,i_c}.n_movie
                AMat = cat(1, AMat, nanmean(res{i,i_c}.interpTracks{i_mov,c},1));
                end
                A{c} = nanmean(AMat,1)+master_adj;
                %if i_c ==1
                    [norm_A_max, it_max] = max(A{c});
                    [norm_A_min, it_min] = min(A{c});
                %end
                
                SEM = nanstd(AMat,[],1)/sqrt(res{i,i_c}.n_movie);
                Amin = A{c} - SEM;
                Aplus = A{c} + SEM;
                    norm_A_max_Amin = Amin(it_max);
                    norm_A_min_Amin = Amin(it_min);
                    norm_A_max_Aplus = Aplus(it_max);
                    norm_A_min_Aplus = Aplus(it_min);
                
                %A{c} = normalize(A{c},'range',[norm_A_min norm_A_max]);
                %Amin = normalize(Amin,'range',[norm_A_min norm_A_max]);
                %Aplus = normalize(Aplus,'range',[norm_A_min norm_A_max]);
                
                iSF = max(AMat);
                res{i,i_c}.A_mean_t{c} = [cT{i_c}{c}-tOffset(c);A{c}];
                res{i,i_c}.A_mean_err{c} = SEM;
                figure(fig_exist{i}); hold on;
                
                if pm.coh_buff
                    ttem = cT{i_c}{c}-tOffset(c);
                Amintem = Amin;
                Aplustem = Aplus;
                Atem = A{c};
                else
                    ttem = cT{i_c}{c}(b+1:end-b)-tOffset(c);
                Amintem = Amin(b+1:end-b);
                Aplustem = Aplus(b+1:end-b);
                Atem = A{c}(b+1:end-b);
                end
                if isempty(pm.ch_col)
                    col_tem=cmap{i_c};
                else
                    col_tem=pm.ch_col(i_ch,:);
                end
                if pm.cohort_norm                  
                    fill([ttem ttem(end:-1:1)], [Amintem Aplustem(end:-1:1)], col_tem,...
                        'EdgeColor', col_tem, 'EdgeAlpha', 0,'FaceAlpha', 0.3,'Parent', ha{i,1});
                    plot(ha{i,1}, ttem, Atem, ip.Results.LineStyle, 'Color', col_tem,...
                        'LineWidth', 1);
                else
                    fill([ttem ttem(end:-1:1)], [Amintem Aplustem(end:-1:1)], col_tem,...
                        'EdgeColor', col_tem, 'EdgeAlpha', 0,'FaceAlpha', 0.3,'Parent', ha{i});
                    plot(ha{i}, ttem, Atem,  ip.Results.LineStyle, 'Color', col_tem,...
                        'LineWidth', 1);
                end
        end
        
        
       %-------------------------------------------------------------------
       if (nCh > 1) 
           i_ch = 2;
       for c = 1:nc{i_c}
                AMat = [];
                for i_mov = 1: res_c{i,i_c}.n_movie
                AMat = cat(1, AMat, nanmean(res_c{i,i_c}.interpTracks{i_mov,c},1));
                end
                A{c} = nanmean(AMat,1);
                SEM = nanstd(AMat,[],1)/sqrt(res_c{i,i_c}.n_movie);
                Amin = A{c} - SEM;
                Aplus = A{c} + SEM;

                
                %norm_A_max = norm_A_max/A{c}(it_max)* max(A{c}) ;
%                 if pm.cohort_norm == true
%                 A{c} = normalize(A{c},'range',[norm_A_min norm_A_max]);
%                 Amin = normalize(Amin,'range',[norm_A_min_Amin norm_A_max_Amin]);
%                 Aplus = normalize(Aplus,'range',[norm_A_min_Aplus norm_A_max_Aplus]);
%                 end

                res_c{i,i_c}.A_mean_t{c} = [cT{i_c}{c}-tOffset(c);A{c}];
                res_c{i,i_c}.A_mean_err{c} = SEM;
                figure(fig_exist{i}); hold on;
                if pm.coh_buff
                    ttem = cT{i_c}{c}-tOffset(c);
                Amintem = Amin;
                Aplustem = Aplus;
                Atem = A{c};
                else
                    ttem = cT{i_c}{c}(b+1:end-b)-tOffset(c);
                Amintem = Amin(b+1:end-b);
                Aplustem = Aplus(b+1:end-b);
                Atem = A{c}(b+1:end-b);
                end

                if isempty(pm.ch_col)
                    col_tem=cmap{i_c};
                else
                    col_tem=pm.ch_col(i_ch,:);
                end
                if pm.cohort_norm                  
                    fill([ttem ttem(end:-1:1)], [Amintem Aplustem(end:-1:1)], col_tem,...
                        'EdgeColor', col_tem, 'EdgeAlpha', 0,'FaceAlpha', 0.3,'Parent', ha{i,2});
                    plot(ha{i,2}, ttem, Atem, '--', 'Color', col_tem,...
                        'LineWidth', 1);
                else
                    fill([ttem ttem(end:-1:1)], [Amintem Aplustem(end:-1:1)], col_tem,...
                        'EdgeColor', col_tem, 'EdgeAlpha', 0,'FaceAlpha', 0.3,'Parent', ha{i});
                    plot(ha{i}, ttem, Atem, '--', 'Color', col_tem,...
                        'LineWidth', 1);
                end
                
       end
       end
       %-------------------------------------------------------------------

    end
end
% for i = 1:num_condition
%     if pm.plot_EpiTIRF == false
%     if pm.cohort_diff == false
%         figure(fig_exist{i}); 
%     else
%         figure(fig_exist{i}); 
%     end
%     else
%     if pm.cohort_diff == false
%         figure(fig_exist{i}); 
%     else
%         figure(fig_exist{i}); 
%     end
%     end
% end
if pm.cohort_diff == false
save([dir_alt filesep 'DAS_cohort.mat'],'res', 'res_c','-v7.3');
end
