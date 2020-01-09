function [fig_exist] = dasPdfConf(data, dir_alt, pm, con_name, DAS_stat, DAS_all, idx, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @iscell);
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('DAS_stat', @iscell);
ip.addRequired('DAS_all', @iscell);
ip.addRequired('idx', @iscell);
ip.addParameter('LT_not_Imax', true, @islogical);
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
ip.parse(data, dir_alt, pm, con_name,DAS_stat, DAS_all, idx, varargin{:});
%==========================================================================
data = ip.Results.data;
pm = ip.Results.pm;
dir_alt=ip.Results.dir_alt;
con_name = ip.Results.con_name;
DAS_stat = ip.Results.DAS_stat;
idx = ip.Results.idx;
DAS_all = ip.Results.DAS_all;
LT_not_Imax = ip.Results.LT_not_Imax;
%==========================================================================
num_condition = max(size(data));
num_movie=zeros(num_condition,1);
for i_condition=1:num_condition
    num_movie(i_condition) = max(size(data{i_condition}));
end
%==========================================================================
% fig_exist = cell(num_condition, 1);
% for i_condition=1:num_condition
%     fig_exist{i_condition} = figure;
% end
fig_exist = cell(1, 2);
fig_exist{1} = figure;
fig_exist{2} = figure;
%==========================================================================
t0=1;
LT = cell(num_condition, pm.num_clus);
Imax = cell(num_condition, pm.num_clus);
for i_condition=1:num_condition
    d_i_c = 1;
    if (pm.plot_visitor == false)
        d_i_c = 2;
    end
    for i_c = 1:d_i_c:pm.num_clus 
        LT{i_condition,i_c} = cell(1,num_movie(i_condition));
        Imax{i_condition,i_c} = cell(1,num_movie(i_condition));     
        for i_mov = 1:num_movie(i_condition)
            id_tem = (idx{i_condition} == i_c) & (DAS_all{i_condition}.MovieNum == i_mov);
            LT{i_condition,i_c}{i_mov} = DAS_all{i_condition}.LT(id_tem);
            Imax{i_condition,i_c}{i_mov} = DAS_all{i_condition}.MaxI(id_tem);
        end
    end
end

curvs = cell(num_condition,1);
for i_condition = 1:num_condition
    if LT_not_Imax
        [fig_pdf, fig_cdf] = plotPDFconf(LT{i_condition,1},LT{i_condition,3});
    else
        [fig_pdf, fig_cdf] = plotPDFconf(Imax{i_condition,1},Imax{i_condition,3},'BoundedSupport', []);
    end
%-------------------------------------------------------------------------------
axes_to_be_copied = findobj(fig_pdf{7},'type','axes');
chilred_to_be_copied = get(axes_to_be_copied,'children');

h = findobj(chilred_to_be_copied,'Type','line');
curvs{i_condition} = [h(1).XData;h(1).YData;h(2).YData;h(3).YData;h(4).XData;h(4).YData;h(5).YData;h(6).YData;];


%-------------------------------------------------------------------------------
close(fig_cdf{:})
close(fig_pdf{:})
end
figure(fig_exist{1}); hold on;
lgd_tem = cell(num_condition*2,1);
for i_condition = 1:num_condition
    figure(fig_exist{1}); hold on;
    plot(curvs{i_condition}(5,:),curvs{i_condition}(6,:),'Color',pm.col_cond(i_condition,:),'Linewidth',1.5);
    figure(fig_exist{2}); hold on;
    plot(curvs{i_condition}(1,:),curvs{i_condition}(2,:),'LineStyle','--','Color',pm.col_cond(i_condition,:),'Linewidth',1.5);
    lgd_tem{i_condition*2-1} = [con_name{i_condition}, ' CCPs'];
    lgd_tem{i_condition*2} = [con_name{i_condition}, ' PCs'];
end

for i_condition = 1:num_condition
    figure(fig_exist{2}); hold on;
    x = curvs{i_condition}(1,:);
    y1 = curvs{i_condition}(3,:);
    y2 = curvs{i_condition}(4,:);
    fill([x fliplr(x)],[y1 fliplr(y2)],pm.col_cond(i_condition,:),'FaceAlpha',0.4,'EdgeAlpha',0);
    plot(curvs{i_condition}(1,:),curvs{i_condition}(2,:),'LineStyle','--','Color',pm.col_cond(i_condition,:),'Linewidth',1.5);
    
    figure(fig_exist{1}); hold on;
    x = curvs{i_condition}(5,:);
    y1 = curvs{i_condition}(7,:);
    y2 = curvs{i_condition}(8,:);
    fill([x fliplr(x)],[y1 fliplr(y2)],pm.col_cond(i_condition,:),'FaceAlpha',0.4,'EdgeAlpha',0);
    plot(curvs{i_condition}(5,:),curvs{i_condition}(6,:),'Color',pm.col_cond(i_condition,:),'Linewidth',1.5);
end
hold off;
%legend(lgd_tem,'Location','northwest');
for i_fig = 1:max(size(fig_exist))
    figure(fig_exist{i_fig});
    if LT_not_Imax
        xlim([5 100]);
        xlabel('$\tau (s)$','Interpreter','latex');
        ylabel('Prob. density');
        set(gca,'Linewidth',1) ;set(gca,'FontSize',8);
    else
        xlim([0 200]);
        xlabel('$I_{max} (a.u.)$','Interpreter','latex');
        ylabel('Prob. density');
        set(gca,'Linewidth',1) ;set(gca,'FontSize',8);
    end
end
