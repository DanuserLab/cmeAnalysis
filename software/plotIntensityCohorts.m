%[cohorts res] = plotIntensityCohorts(data, varargin) displays average CCP intensities for a range of lifetime cohorts
%
% Inputs:
%         data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%  'MaxIntensityThreshold' : threshold value generated by runLifetimeAnalysis()
%              'SlaveName' : cell array of strings describing slave channels, for legends
%          'ScalingFactor' : scaling vector to adjust between channel intensities
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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

% Francois Aguet (last modified 04/30/2013)

function [cohorts, res, ha] = plotIntensityCohorts(data, varargin)

nCh = numel(data(1).channels);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', nCh:-1:1);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('ShowVariation', true, @islogical);
ip.addParamValue('FillMode', 'SEM', @(x) any(strcmpi(x, {'SEM', 'pct'})));
ip.addParamValue('FrontLayer', false, @islogical);
ip.addParamValue('ShowBackground', false, @islogical);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('RescalingReference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
ip.addParamValue('ScaleSlaveChannel', true, @islogical);
ip.addParamValue('ScalingFactor', [], @(x) numel(x)==nCh);
ip.addParamValue('MaxIntensityThreshold', 0);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('SlaveName', [], @iscell);
ip.addParamValue('ChannelNames', []);
ip.addParamValue('LineStyle', '-');
ip.addParamValue('Hues', []);
ip.addParamValue('Colormap', []);
ip.addParamValue('ColormapFill', []);
ip.addParamValue('DisplayMode', 'screen');
ip.addParamValue('DisplayAll', false);
ip.addParamValue('TrackIndex', []);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Alpha', 0.05);
ip.addParamValue('YTick', []);
ip.addParamValue('YLim', []);
ip.addParamValue('Parent', []);
ip.addParamValue('RemoveOutliers', false, @islogical);
ip.addParamValue('PlotXLabel', true, @islogical);
ip.addParamValue('ShowPct', true, @islogical);
ip.addParamValue('ShowStats', false, @islogical);
ip.addParamValue('AvgFun', @nanmean, @(x) isa(x, 'function_handle'));
ip.addParamValue('LftDataName', 'lifetimeData.mat');
ip.addParamValue('AmplitudeCorrection', []);
ip.addParamValue('Align', 'left', @(x) any(strcmpi(x, {'right', 'left'})));
% ip.addParamValue('MinTracksPerCohort', 5);
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;
sf = ip.Results.ScalingFactor;
hues = ip.Results.Hues;

% load data
lftData = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite,...
    'LifetimeData', ip.Results.LftDataName, 'Scale', ip.Results.Rescale,...
    'Cutoff_f', ip.Results.Cutoff_f, 'ReturnValidOnly', true,...
    'ExcludeVisitors', ip.Results.ExcludeVisitors, 'Mask', true,...
    'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);

% if no specific channel is selected, all channels are shown
chVec = ip.Results.ch;
mCh = find(strcmp(data(1).source, data(1).channels));

nd = numel(data);
kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

nc = numel(cohortBounds)-1;
b = size(lftData(1).sbA,2);
framerate = data(1).framerate;

% # data points in cohort (including buffer frames)
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);

% time vectors for cohorts
if strcmpi(ip.Results.Align, 'left')
    cT = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'unif', 0);
else
    cT = arrayfun(@(i) (-i+b+1:b)*framerate, iLength, 'unif', 0);
end
XLim = [cT{end}(1)-b cT{end}(end)+b];

YLim = ip.Results.YLim;
if isempty(YLim) && ~isempty(ip.Results.YTick)
    YLim = ip.Results.YTick([1 end]);
end


if ~isempty(ip.Results.TrackIndex)
    lftFields = fieldnames(lftData);
    i = cellfun(@(f) size(lftData(1).(f),1)==numel(lftData(1).lifetime_s), lftFields);
    lftFields = lftFields(i);
    for i = 1:nd
        for f = 1:numel(lftFields)
            lftData(i).(lftFields{f}) = lftData(i).(lftFields{f})(ip.Results.TrackIndex{i},:,:);
        end
    end
end

% loop through data sets, generate cohorts for each
res(1:nd) = struct('interpTracks', [], 'interpSigLevel', []);
for i = 1:nd
    
    % for intensity threshold in master channel
    maxA = max(lftData(i).A(:,:,mCh), [], 2);
    
    for ch = 1:nCh % channels
        % interpolate tracks to mean cohort length
        for c = 1:nc % cohorts
            % tracks in current cohort (above threshold)
            if c<nc
                cidx = find(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1) &...
                    maxA > ip.Results.MaxIntensityThreshold);
            else % inclusive upper bound for last cohort
                cidx = find(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<=cohortBounds(c+1) &...
                    maxA > ip.Results.MaxIntensityThreshold);
            end
            
            nt = numel(cidx);
            if nt>0
                interpTracks = zeros(nt,iLength(c));
                sigma_rMat = zeros(nt,iLength(c));
                cLengths = lftData(i).trackLengths(cidx);
                % loop through track lengths within cohort
                for t = 1:nt
                    A = [lftData(i).sbA(cidx(t),:,ch) lftData(i).A(cidx(t),1:cLengths(t),ch) lftData(i).ebA(cidx(t),:,ch)];
                    bgr = [lftData(i).sbSigma_r(cidx(t),:,ch) lftData(i).sigma_r(cidx(t),1:cLengths(t),ch) lftData(i).ebSigma_r(cidx(t),:,ch)];
                    A(isnan(A)) = nanmin(A);
                    bgr(isnan(bgr)) = nanmin(bgr);
                    % align to track start
                    %w = min(numel(A),iLength);
                    %interpTracks(t,1:w) = A(1:w);
                    %sigma_r_Ia(t,1:w) = bgr(1:w);
                    
                    % interpolate to mean length
                    xi = linspace(1,cLengths(t)+2*b, iLength(c));
                    interpTracks(t,:) = binterp(A, xi);
                    sigma_rMat(t,:) = binterp(bgr, xi);
                    %interpTracks(t,:) = interp1(1:cLengths(t)+2*b, A, xi, 'cubic');
                    %sigma_rMat(t,:) = interp1(1:cLengths(t)+2*b, bgr, xi, 'cubic');
                end
                
                res(i).interpTracks{ch,c} = interpTracks;
                res(i).interpSigLevel{ch,c} = kLevel*sigma_rMat;
                % split as a function of slave channel signal
                if isfield(lftData(i), 'significantMaster') && nCh>1
                    sigIdx = lftData(i).significantMaster(:,ch)==1;
                    %sigIdx = lftData(i).significantSlave(:,ch)==1;
                    %sigIdx = lftData(i).significantMaster(:,ch)==0 & lftData(i).significantSlave(:,ch)==1;
                    res(i).sigIdx{c}(:,ch) = sigIdx(cidx);
                else
                    res(i).sigIdx{c}(:,ch) = ones(numel(cidx),1);
                end
                if ch==1
                    res(i).trackIdx{c} = lftData(i).index(cidx);
                end
            else
                res(i).interpTracks{ch,c} = NaN(1,iLength(c));
                res(i).interpSigLevel{ch,c} = NaN(1,iLength(c));
                res(i).sigIdx{c}(:,ch) = NaN;
                if ch==1
                    res(i).trackIdx{c} = NaN;
                end
            end
        end
    end
end

cohortLabels = arrayfun(@(i) [num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)) ' s'], 1:nc, 'Unif', 0);
XTick = (cohortBounds(1:end-1)+[cohortBounds(2:end-1) cohortBounds(end)-framerate])/2;
if strcmpi(ip.Results.Align, 'right')
    XTick = 0 - XTick(end:-1:1);
    cohortLabels = cohortLabels(end:-1:1);
end


fset = loadFigureSettings(ip.Results.DisplayMode);

% Set colormap depending on # channels
cmap = ip.Results.Colormap;
if isempty(cmap)
    cmap = cell(1,nCh);
    if nCh==1
        if isempty(hues)
            hues = getFluorophoreHues(data(1).markers);
        end
        v = mod(hues(1)+linspace(-0.05, 0.05, nc)', 1);
        cmap{1} = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
    else
        if isempty(hues)
            hues = getFluorophoreHues(data(1).markers);
        end
        if nCh==2
            hb = 0.1;
        else
            hb = 0.05;
        end
        for ch = 1:nCh
            v = mod(hues(ch)+linspace(-hb, hb, nc)', 1);
            cmap{ch} = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
        end
    end
end
cv = ip.Results.ColormapFill;
if isempty(cv)
    cv = cell(1,nCh);
    for c = 1:nCh
        tmp = rgb2hsv(cmap{c});
        cv{c} = hsv2rgb([tmp(:,1) 0.4*ones(nc,1) ones(nc,1)]);
    end
end
% scale slave channels relative to master (for visualization only)
if isempty(sf) 
    if ip.Results.ScaleSlaveChannel && nCh>1
        for ch = 1:nCh
            iSF = zeros(1,nc);
            for c = 1:nc
                % find largest mean of all cohorts
                M = arrayfun(@(x) mean(x.interpTracks{ch,c},1), res, 'unif', 0);
                M = mean(vertcat(M{:}), 1);
                iSF(c) = max(M);
            end
            sf(ch) = max(iSF);
        end
        sf = sf(mCh)./sf;
    else
        sf = ones(1,nCh);
    end
end        

% output
cohorts.bounds = cohortBounds;



%==================================================
% Plot cohorts
%==================================================
switch nCh
    case 1
        ah = 1;
        na = 1;
        sigCombIdx = [];
    case 2
        na = 2;
        ah = 1;
        sigCombIdx = [1 0]';
        
        pct = zeros(nd,2);
        for i = 1:nd
            vidx = max(lftData(i).A(:,:,mCh),[],2) > ip.Results.MaxIntensityThreshold;
            s = lftData(i).significantMaster(vidx,:);
            %pct(i,:) = sum([s(:,2) ~s(:,2)],1)/size(s,1);
            idx = lftData(i).maxA(vidx,1)>ip.Results.MaxIntensityThreshold;
            pct(i,:) = sum([s(idx,2) ~s(idx,2)],1)/sum(idx);
        end
        meanPct = mean(pct,1);
        stdPct = std(pct,[],1);
    case 3
        na = 4;
        ah = 2;
        sigCombIdx = [1 1; 1 0; 0 1; 0 0];
        
        pct = zeros(nd,4);
        for i = 1:nd
            vidx = max(lftData(i).A(:,:,mCh),[],2) > ip.Results.MaxIntensityThreshold;
            s = lftData(i).significantMaster(vidx,:);
            pct(i,:) = sum([s(:,2)&s(:,3) s(:,2)&~s(:,3) ~s(:,2)&s(:,3) ~s(:,2)&~s(:,3)],1)/size(s,1);
        end
        meanPct = mean(pct,1);
        stdPct = std(pct,[],1);
end

if ~isempty(sigCombIdx)
    SlaveName = ip.Results.SlaveName;
    if isempty(SlaveName)
        SlaveName = data(1).markers(2:nCh);
    end
    tmp = sigCombIdx;
    tmp(tmp==1) = '+';
    tmp(tmp==0) = '-';
    atext = cell(1,na);
    switch nCh
        case 2
            for a = 1:na
                atext{a} = [tmp(a,1) SlaveName{1} ': ' num2str(meanPct(a)*100, '%.1f') '�' num2str(stdPct(a)*100, '%.1f') '%'];
            end
        case 3
            for a = 1:na
                atext{a} = [SlaveName{1} tmp(a,1) ' / ' SlaveName{2} tmp(a,2) ': '...
                    num2str(meanPct(a)*100, '%.1f') '�' num2str(stdPct(a)*100, '%.1f') '%'];
            end
    end
end

aw = ceil(na/ah);
aposy = 2;


A = cell(nCh,nc);
ha = ip.Results.Parent;
if isempty(ha)
    if ip.Results.ShowPct && nCh>2
        ha = setupFigure(ah, aw, 'YSpace', [3 1 0.5], 'XSpace', [2 0.5 3.5],...
            'SameAxes', true, 'Name', 'Intensity cohorts', 'DisplayMode', ip.Results.DisplayMode);
    else
        ha = setupFigure(ah, aw, 'YSpace', [2.5 1 1], 'XSpace', [2 0.75 0.5],...
            'SameAxes', true, 'Name', 'Intensity cohorts', 'DisplayMode', ip.Results.DisplayMode);
    end
end

% now plot cohorts for each combination
if strcmpi(ip.Results.Align, 'left')
    tOffset = @(c) 0;
else
    tOffset = @(c) cT{c}(end)-b*framerate;
end
for a = 1:na

    % combination for these axes: sigCombIdx(a)
    for i = 1:nd
        for c = 1:nc
            switch nCh
                case 1
                    res(i).sigComb{a,c} = res(i).sigIdx{c}==1;
                case 2
                    res(i).sigComb{a,c} = res(i).sigIdx{c}(:,2)==sigCombIdx(a,1);
                case 3
                    res(i).sigComb{a,c} = res(i).sigIdx{c}(:,2)==sigCombIdx(a,1) &...
                        res(i).sigIdx{c}(:,3)==sigCombIdx(a,2);
            end
        end
    end
    
    for c = nc:-1:1
        for ch = chVec; % plot master channel last
            if nd > 1
                % means for each data set
                AMat = arrayfun(@(x) ip.Results.AvgFun(x.interpTracks{ch,c}(x.sigComb{a,c},:),1), res, 'unif', 0);
                AMat = vertcat(AMat{:});
                % # of tracks from each data set in this cohort
                ntCoSel = arrayfun(@(x) sum(x.sigComb{a,c}), res);
                A{ch,c} = nanmean(AMat,1);
                SEM = nanstd(AMat,[],1)/sqrt(nd);
                Amin = A{ch,c} - SEM;
                Aplus = A{ch,c} + SEM;
            else
                % if input is a single data set, show median + percentiles
                M = prctile(res(1).interpTracks{ch,c}(res(1).sigComb{a,c},:), [25 50 75], 1);
                ntCoSel = sum(res(1).sigComb{a,c});
                A{ch,c} = M(2,:);
                Amin = M(1,:);
                Aplus = M(3,:);
            end
            % plot cohort only if at least half the data sets have tracks in this cohort
            if sum(ntCoSel>0)/numel(ntCoSel) > 0.5
                if ip.Results.ShowVariation
                    fill([cT{c} cT{c}(end:-1:1)]-tOffset(c), sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:),...
                        'EdgeColor', cmap{ch}(c,:), 'Parent', ha(a));
                end
                if ~ip.Results.FrontLayer
                    plot(ha(a), cT{c}-tOffset(c), sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:),...
                        'LineWidth', 1);
                end
            end
            cohorts(a).t{c} = cT{c};
            cohorts(a).Amin{ch,c} = Amin;
            cohorts(a).Aplus{ch,c} = Aplus;
            cohorts(a).A{ch,c} = A{ch,c};
        end
    end
    
    for ch = chVec
        % Plot mean/median in front
        if ip.Results.FrontLayer
            for c = nc:-1:1
                plot(ha(a), cT{c}-tOffset(c), sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
            end
        end
        
        % Plot signifcance threshold in front
        if ip.Results.ShowBackground && ch~=mCh
            % Background level: median of all detections
            if nd>1
                % median background level per cohort for each data set
                medM = arrayfun(@(i) cellfun(@(x) nanmean(x(:)), i.interpSigLevel(ch,:)) , res, 'unif', 0);
                medM = vertcat(medM{:});
                plot(ha(a), [-10 120], sf(ch)*nanmean(medM(:))*[1 1], 'k--', 'LineWidth', 1);
            else
                % median background level per cohort
                medC = cellfun(@(x) nanmean(x(:)), res.interpSigLevel(ch,:));
                plot(ha(a), [-10 120], nanmean(medC)*[1 1], 'k--', 'LineWidth', 1);
            end
        end
    end
end

if isempty(YLim)
    YLim = get(ha, 'YLim');
    if na>1
        YLim = vertcat(YLim{:});
        YLim = [min(YLim(:,1)) max(YLim(:,2))];
    end
end

if ~isempty(sigCombIdx)
    for a = 1:na
        text(XLim(2), YLim(2), atext{a}, fset.sfont{:},...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Parent', ha(a));
    end
end

if ip.Results.ShowPct && nCh>2
    dy = 0.75;
    hav = axes(fset.axOpts{:}, 'Position', [1.5+ceil(na/ah)*6.75 aposy+ah*3.5+(ah-1)*dy-1.6 2.4 1.6], 'Units', 'normalized');
    vennplot(meanPct(2), meanPct(3), meanPct(1), ip.Results.SlaveName,...
        'Parent', hav, 'Font', fset.sfont);
    axis off;
end
        
%     if ip.Results.ShowPct
%         axes(fset.axOpts{:}, 'Position', [15.5 2 3 2.5], 'TickLength', fset.TickLength*6/3);
%         barplot2(mean(M,1)', std(M,[],1)', 'Angle', 0, 'BarWidth', 1, 'GroupDistance', 1,...
%             'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.4*[1 1 1],...
%             'YLim', [0 100], 'LineWidth', 1);
%         set(gca, 'FontSize', 8);
%         
%         h = title(['% ' ip.Results.SlaveName ' pos. CCPs'], fset.sfont{:});
%         %h = ylabel('% CCPs/cohort', fset.lfont{:});
%         pos = get(h, 'Position');
%         %pos(1) = 0.8*pos(1);
%         pos(2) = 1.1*pos(2);
%         set(h, 'Position', pos);
%         set(gca, 'YTick', 0:20:100, 'XTickLabel', cohortLabels);
%         rotateXTickLabels(gca, 'AdjustFigure', false);
%         xlabel('Lifetime cohort', fset.lfont{:});
%     end
    idx = 1-mod(1:na,2)~=0;
    set(ha(idx), 'YTickLabel', []);
    arrayfun(@(x) ylabel(x, 'Fluo. intensity (A.U.)', fset.lfont{:}), ha(~idx));
% end

set(ha, 'XLim', XLim, 'XTick', XTick, 'YLim', YLim);
if ~isempty(ip.Results.YTick)
    set(ha, 'YTick', ip.Results.YTick);
end

if ip.Results.PlotXLabel
    idx = ah-ceil((1:na)/2)>0;
    set(ha, 'XTick', XTick, 'XTickLabel', cohortLabels);
    set(ha(idx), 'XTickLabel', []);
    for i = find(~idx)
        rotateXTickLabels(ha(i), 'AdjustFigure', false);
        xlabel(ha(i), 'Lifetime cohort', fset.lfont{:});
    end
end


if ip.Results.ShowStats

    % plot total tracks in each cohort
    for a = 1:na
        % events selected/cohort
        ntCoSel = arrayfun(@(c) arrayfun(@(x) sum(x.sigComb{a,c}), res), 1:nc, 'unif', 0);
        ntCoSel = vertcat(ntCoSel{:})';
        
        %text(XLim(1)*0.95, 1.0*YLim(end), '# CCPs:', 'Parent',  ha(a), fset.sfont{:},...
        %    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        for c = 1:nc
            %if sum(ntCoSel{c}>0)/numel(ntCoSel{c}) > 0.5
                text(XTick(c), YLim(end), num2str(sum(ntCoSel(:,c))), 'Parent',  ha(a), fset.sfont{:},...
                   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            %end
        end
        
%         % percentages
%         pct = ntCoSel./repmat(sum(ntCoSel,2), [1 nc])*100;
%         text(XLim(1)*0.95, 1.0*YLim(end), '% CCPs:', 'Parent',  ha(a), fset.sfont{:},...
%             'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
%         for c = 1:nc
%                 text(XTick(c), 1.0*YLim(end), num2str(mean(pct(:,c),1), '%.0f'), 'Parent',  ha(a), fset.sfont{:},...
%                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%         end        
    end
end
