

function [figH] = dasFigure1(data, DAS_all,Track_info, DAS_stat, idx, con_name,dir_alt, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) iscell(x));
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addRequired('DAS_stat', @(x) iscell(x));
ip.addRequired('idx', @(x) iscell(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('fix_cluster_num', 2, @isnumeric);
ip.addParameter('i_condition', 1, @isnumeric);
ip.addParameter('col_scatter', [0 1 1], @isnumeric);
ip.addParameter('fig_name', 'fig1', @ischar);
ip.addParameter('col_scatter_edge', [0 0 1], @isnumeric);
ip.parse(data, DAS_all, Track_info,DAS_stat, idx, con_name,dir_alt, pm, varargin{:});


DAS_all = ip.Results.DAS_all;
con_name = ip.Results.con_name;
num_condition = max(size(DAS_all));
idx = ip.Results.idx;
data = ip.Results.data;
dir_alt = ip.Results.dir_alt;
Track_info = ip.Results.Track_info;
DAS_stat = ip.Results.DAS_stat;
i_condition = ip.Results.i_condition;
fig_name = ip.Results.fig_name;
pm = ip.Results.pm;
col_scatter = ip.Results.col_scatter;
col_scatter_edge = ip.Results.col_scatter_edge;
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
figH = cell(3,1);
%==========================================================================
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.addGroup(fig_name);
desktop.setGroupDocked(fig_name, 0);
myDim   = java.awt.Dimension(3, 3);   % 4 columns, 2 rows
% 1: Maximized, 2: Tiled, 3: Floating
n_fig = 9;
figH{1}    = gobjects(n_fig, 1);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
for iFig = 1:n_fig
   figH{1}(iFig) = figure('WindowStyle', 'docked', ...
      'Name', sprintf('Figure %d', iFig), 'NumberTitle', 'off');
   set(get(handle(figH{1}(iFig)), 'javaframe'), 'GroupName', fig_name);
end
warning(bakWarn);
desktop.setDocumentArrangement(fig_name, 2, myDim)
%==========================================================================
num_clus =3;
%==========================================================================
label1 = pm.label1;
label2 = pm.label2;
label3 = pm.label3;
%==========================================================================
color_clus = pm.color_clus;
%==========================================================================
%==========================================================================
for iFig = 1:n_fig
   figure(figH{1}(iFig));
end
%==========================================================================
figure(figH{1}(2));
if pm.fig1_LT_Imax_stack == false
line(DAS_stat{i_condition}.Imax,DAS_stat{i_condition}.pdf_Imax,'Color','r','Linewidth',1)
xlim([0 200])
ylabel('Prob. Density'); xlabel('$I_{max}$ (a.u.)','interpreter','latex')
set(gca,'Linewidth',2) 
set(gca,'FontSize',8)
%set(gca,'fontsize',fontsize_def, 'FontWeight', 'bold')

ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
set(ax1,'Position',ax1.Position+[0 0 -0.02 0])
ax1_pos = ax1.Position; % position of first axes

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
set(gca,'Linewidth',1) 
set(gca,'FontSize',8)

line(DAS_stat{i_condition}.LT,DAS_stat{i_condition}.pdf_LT,'Parent',ax2,'Color','k','Linewidth',1);
xlim([5 150])
set(gca,'Linewidth',1) 
ylabel('Prob. Density');xlabel('$\tau$ (s)','interpreter','latex')
else
    subplot(2,1,2);
    line(DAS_stat{i_condition}.Imax,DAS_stat{i_condition}.pdf_Imax,'Color','k','Linewidth',1);
    xlim([0 200]);
    ylabel('Prob. Density'); xlabel('$I_{max}$ (a.u.)','interpreter','latex');
    set(gca,'Linewidth',1); set(gca,'FontSize',8);
    subplot(2,1,1);
    line(DAS_stat{i_condition}.LT,DAS_stat{i_condition}.pdf_LT,'Color','k','Linewidth',1);
    xlim([5 150])
    set(gca,'Linewidth',1); set(gca,'FontSize',8);
    ylabel('Prob. Density');xlabel('$\tau$ (s)','interpreter','latex')
end
%==========================================================================
figure(figH{1}(4));
line(DAS_stat{i_condition}.d,DAS_stat{i_condition}.pdf_d,'Color','k','Linewidth',1)
ylabel('Prob. Density');xlabel(label2,'interpreter','latex')
set(gca,'Linewidth',1) 
xlim([-1.5 1])
figure(figH{1}(5));
line(DAS_stat{i_condition}.dv,DAS_stat{i_condition}.pdf_dv,'Color','k','Linewidth',1)
ylabel('Prob. Density');xlabel(label1,'interpreter','latex')
set(gca,'Linewidth',1) 
xlim([-5 0])
figure(figH{1}(6));
line(DAS_stat{i_condition}.d32,DAS_stat{i_condition}.pdf_d32,'Color','k','Linewidth',1)
ylabel('Prob. Density');xlabel(label3,'interpreter','latex')
set(gca,'Linewidth',1) 
xlim([-3 3])
%==========================================================================
%==========================================================================
dasClusterPlot_single(DAS_stat{i_condition}.z_d_dv, DAS_stat{i_condition}.edges_d_dv, {label2,label1},'fig_exist',figH{1}(7));
dasClusterPlot_single(DAS_stat{i_condition}.z_d_d32, DAS_stat{i_condition}.edges_d_d32, ...
    {label1,label3},'fig_exist',figH{1}(8));
dasClusterPlot_single(DAS_stat{i_condition}.z_d32_dv, DAS_stat{i_condition}.edges_d32_dv, ...
    {label3,label2},'fig_exist',figH{1}(9));
%==========================================================================
% dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_var, DAS_all{i_condition}.DAS , idx{i_condition},DAS_stat{i_condition}.edges_d_dv, {label1,label2},'fig_exist',figH{1}(10));
% dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3, DAS_all{i_condition}.DAS , idx{i_condition}, DAS_stat{i_condition}.edges_d_d32,...
%     {label1,label3},'fig_exist',figH{1}(11));
% dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_var, DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3 , idx{i_condition}, DAS_stat{i_condition}.edges_d32_dv,...
%     {label3,label2},'fig_exist',figH{1}(12));
%==========================================================================
%==========================================================================
figure(figH{1}(3));

dasComputeD(data{1}, Track_info,[1 0],pm,DAS_all{1}.ImaxAll,'plotD', true,'fig_exist',figH{1}(3));
hold on;
id_temp = (DAS_all{1}.MovieNum==1) & (DAS_all{1}.LT==30) & (idx{1} == 1);
id_temp = DAS_all{1}.TrackID(id_temp);
i_tem=randsample(max(size(id_temp)),1);
plot(1:30,Track_info{1}.A(id_temp(i_tem),1:30),'linewidth',1,'color',pm.color_clus{1},...
    'Marker','o',...
    'LineWidth',1,...
    'MarkerSize',1,...
    'MarkerEdgeColor',pm.color_clus{1},...
    'MarkerFaceColor',pm.color_clus{1})
% Track_info{1}.A(id_temp(6),1:lt_tem) for paper (Zhiming 171219 ctrl)
id_temp = (DAS_all{1}.MovieNum==1) & (DAS_all{1}.LT==10) & (idx{1} == 3);
id_temp = DAS_all{1}.TrackID(id_temp);
hold on;
i_tem=randsample(max(size(id_temp)),1);
plot(1:10,Track_info{1}.A(id_temp(i_tem),1:10),'linewidth',1,'color',pm.color_clus{3},...
    'Marker','o',...
    'LineWidth',1,...
    'MarkerSize',1,...
    'MarkerEdgeColor',pm.color_clus{3},...
    'MarkerFaceColor',pm.color_clus{3})
% Track_info{1}.A(id_temp(3),1:lt_tem) for paper (Zhiming 171219 ctrl)
lt_tem = 10;
id_temp = (DAS_all{1}.MovieNum==1) & (DAS_all{1}.LT==lt_tem) & (idx{1} == 2);
id_temp = DAS_all{1}.TrackID(id_temp);
hold on;
i_tem=randsample(max(size(id_temp)),1);
plot(1:lt_tem,Track_info{1}.A(id_temp(i_tem),1:lt_tem),'linewidth',1,'color',pm.color_clus{2},...
    'Marker','o',...
    'LineWidth',1,...
    'MarkerSize',1,...
    'MarkerEdgeColor',pm.color_clus{2},...
    'MarkerFaceColor',pm.color_clus{2})
% Track_info{1}.A(id_temp(8),1:lt_tem) for paper (Zhiming 171219 ctrl)
%==========================================================================
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.addGroup([fig_name '_2']);
desktop.setGroupDocked([fig_name '_2'], 0);
myDim   = java.awt.Dimension(4, 2);   % 4 columns, 2 rows
% 1: Maximized, 2: Tiled, 3: Floating
n_fig = 8;
figH{2}    = gobjects(n_fig, 1);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
for iFig = 1:n_fig
   figH{2}(iFig) = figure('WindowStyle', 'docked', ...
      'Name', sprintf('Figure %d', iFig), 'NumberTitle', 'off');
   set(get(handle(figH{2}(iFig)), 'javaframe'), 'GroupName', [fig_name '_2']);
end
warning(bakWarn);
desktop.setDocumentArrangement([fig_name '_2'], 2, myDim)
%==========================================================================
%==========================================================================
for iFig = 1:n_fig
   figure(figH{2}(iFig));
end
%==========================================================================
Z_norm=[DAS_all{i_condition}.DAS_var,DAS_all{i_condition}.DAS];

q_low_LT = 0.02; q_high_LT = 0.96; q_low_Imax = 0.01; q_high_Imax = 0.90;
Imax_range = zeros(num_clus,2);
LT_range = zeros(num_clus,2);
for i_c = 1: num_clus
    id_temp = (idx{i_condition}==i_c);
    Imax_range(i_c,1) = quantile(DAS_all{i_condition}.MaxI(id_temp),q_low_Imax);
    Imax_range(i_c,2) = quantile(DAS_all{i_condition}.MaxI(id_temp),q_high_Imax);
    LT_range(i_c,1) = quantile(DAS_all{i_condition}.LT(id_temp),q_low_LT);
    LT_range(i_c,2) = quantile(DAS_all{i_condition}.LT(id_temp),q_high_LT);
end

Imax_range(2,1) = 0; Imax_range(2,2) = 1000000; LT_range(2,1) = 0; LT_range(2,2) = 1000000;

        id_Imax = (DAS_all{i_condition}.MaxI > max(Imax_range(:,1))) & (DAS_all{i_condition}.MaxI < min(Imax_range(:,2)));
        id_LT = (DAS_all{i_condition}.LT > max(LT_range(:,1))) & (DAS_all{i_condition}.LT < min(LT_range(:,2)));
        hold on;

%==========================================================================


figure(figH{2}(1));xlabel('$\tau$ (s)','interpreter','latex'); ylabel('Prob. Density');
line(DAS_stat{i_condition}.LT,DAS_stat{i_condition}.pdf_LT,'Color','k','Linewidth',2); hold on
figure(figH{2}(3));xlabel('$I_{max}$ (a.u.)','interpreter','latex'); ylabel('Prob. Density');
line(DAS_stat{i_condition}.Imax,DAS_stat{i_condition}.pdf_Imax,'Color','k','Linewidth',2); hold on

    figure(figH{2}(1));  
    xlim([5 100]);
    line(DAS_stat{i_condition}.LT1,DAS_stat{i_condition}.pdf_LT1*DAS_stat{i_condition}.prop_c(1),'Color',color_clus{1},'Linewidth',2); hold on
    line(DAS_stat{i_condition}.LT2,DAS_stat{i_condition}.pdf_LT2*DAS_stat{i_condition}.prop_c(2),'Color',color_clus{2},'Linewidth',2); hold on
    line(DAS_stat{i_condition}.LT3,DAS_stat{i_condition}.pdf_LT3*DAS_stat{i_condition}.prop_c(3),'Color',color_clus{3},'Linewidth',2); hold on

    
    figure(figH{2}(3));
    xlim([0 200]);
    line(DAS_stat{i_condition}.Imax1,DAS_stat{i_condition}.pdf_Imax1*DAS_stat{i_condition}.prop_c(1),'Color',color_clus{1},'Linewidth',2); hold on
    line(DAS_stat{i_condition}.Imax2,DAS_stat{i_condition}.pdf_Imax2*DAS_stat{i_condition}.prop_c(2),'Color',color_clus{2},'Linewidth',2); hold on
    line(DAS_stat{i_condition}.Imax3,DAS_stat{i_condition}.pdf_Imax3*DAS_stat{i_condition}.prop_c(3),'Color',color_clus{3},'Linewidth',2); hold on
    


%==========================================================================
y_LT_min = zeros(size(DAS_stat{i_condition}.LT));
LT_min1 = min(DAS_stat{i_condition}.LT1); LT_max1 = max(DAS_stat{i_condition}.LT1);
LT_min3 = min(DAS_stat{i_condition}.LT3); LT_max3 = max(DAS_stat{i_condition}.LT3);
for i = 1:max(size(DAS_stat{i_condition}.LT))
    if (i >= LT_min1) && (i >= LT_min3) && (i <= LT_max1) && (i <= LT_max3) 
        [~, i1] = min(abs(DAS_stat{i_condition}.LT1-DAS_stat{i_condition}.LT(i)));
        [~, i3] = min(abs(DAS_stat{i_condition}.LT3-DAS_stat{i_condition}.LT(i)));
        y_LT_min(i) = min([DAS_stat{i_condition}.pdf_LT1(i1)*DAS_stat{i_condition}.prop_c(1),DAS_stat{i_condition}.pdf_LT3(i3)*DAS_stat{i_condition}.prop_c(3)]);
    end
end
%==========================================================================
figure(figH{2}(1));
id_temp = (DAS_stat{i_condition}.LT1 >= max(LT_range(:,1))) & (DAS_stat{i_condition}.LT1 <= min(LT_range(:,2)));
area(DAS_stat{i_condition}.LT1(id_temp),DAS_stat{i_condition}.pdf_LT1(id_temp)*DAS_stat{i_condition}.prop_c(1),'FaceColor', [0.5 0.5 0.5],'FaceAlpha',0.5, 'EdgeColor','none');
id_temp = (DAS_stat{i_condition}.LT3 >= max(LT_range(:,1))) & (DAS_stat{i_condition}.LT3 <= min(LT_range(:,2)));
area(DAS_stat{i_condition}.LT3(id_temp),DAS_stat{i_condition}.pdf_LT3(id_temp)*DAS_stat{i_condition}.prop_c(3),'FaceColor', [0.5 0.5 0.5],'FaceAlpha',0.5, 'EdgeColor','none');
ylim([0 0.025])
%==========================================================================
y_Imax_min = zeros(size(DAS_stat{i_condition}.Imax));
Imax_min1 = min(DAS_stat{i_condition}.Imax1); Imax_max1 = max(DAS_stat{i_condition}.Imax1);
Imax_min3 = min(DAS_stat{i_condition}.Imax3); Imax_max3 = max(DAS_stat{i_condition}.Imax3);
for i = 1:max(size(DAS_stat{i_condition}.Imax))
    if (i >= Imax_min1) && (i >= Imax_min3) && (i <= Imax_max1) && (i <= Imax_max3) 
        [~, i1] = min(abs(DAS_stat{i_condition}.Imax1-DAS_stat{i_condition}.Imax(i)));
        [~, i3] = min(abs(DAS_stat{i_condition}.Imax3-DAS_stat{i_condition}.Imax(i)));
        y_Imax_min(i) = min([DAS_stat{i_condition}.pdf_Imax1(i1)*DAS_stat{i_condition}.prop_c(1),DAS_stat{i_condition}.pdf_Imax3(i3)*DAS_stat{i_condition}.prop_c(3)]);
    end
end
%==========================================================================
figure(figH{2}(3));
id_temp = (DAS_stat{i_condition}.Imax1 >= max(Imax_range(:,1))) & (DAS_stat{i_condition}.Imax1 <= min(Imax_range(:,2)));
area(DAS_stat{i_condition}.Imax1(id_temp),DAS_stat{i_condition}.pdf_Imax1(id_temp)*DAS_stat{i_condition}.prop_c(1),'FaceColor', [0.5 0.5 0.5],'FaceAlpha',0.5, 'EdgeColor','none');
id_temp = (DAS_stat{i_condition}.Imax3 >= max(Imax_range(:,1))) & (DAS_stat{i_condition}.Imax3 <= min(Imax_range(:,2)));
area(DAS_stat{i_condition}.Imax3(id_temp),DAS_stat{i_condition}.pdf_Imax3(id_temp)*DAS_stat{i_condition}.prop_c(3),'FaceColor', [0.5 0.5 0.5],'FaceAlpha',0.5, 'EdgeColor','none');
ylim([0 0.02])
%==========================================================================
figure(figH{2}(5));
num_track_select = 20;
size_exp = 5;
map_col = 'summer';

dist_temp = sqrt((Z_norm(:,2)-DAS_stat{i_condition}.d_mod(1)).^2+(Z_norm(:,1)-DAS_stat{i_condition}.dv_mod(1)).^2);
id_temp = id_LT & (DAS_all{i_condition}.MovieNum==1) & (idx{i_condition} == 1);
DAS_track_id = DAS_all{i_condition}.TrackID(id_temp);
DAS_dist = dist_temp(id_temp);
Z1=Z_norm(id_temp,1);
Z2=Z_norm(id_temp,2);
[~, id_sort] = sort(DAS_dist);
id = DAS_track_id(id_sort);
id = id(1:num_track_select);
Z1 = Z1(id_sort);Z1=Z1(1:num_track_select);
Z2 = Z2(id_sort);Z2=Z2(1:num_track_select);
plotExampleTrack(data{i_condition}, Track_info,id,{gcf},pm,'subplot_id',[1,2,1:2]);
if pm.fig1_mosaic
subplot(1,2,2); 
end
xlim([0 30]); ylim([0 120]);
figure(figH{2}(2));
dasClusterPlot_single(DAS_stat{i_condition}.z_d_dv, DAS_stat{i_condition}.edges_d_dv, {label2,label1},'fig_exist',figH{2}(2));
hold on
scatter(Z2,Z1,'filled','SizeData',size_exp,'MarkerFaceColor',col_scatter,'MarkerEdgeColor',col_scatter_edge);alpha(1);

figure(figH{2}(6));
dist_temp = sqrt((Z_norm(:,2)-DAS_stat{i_condition}.d_mod(3)).^2+(Z_norm(:,1)-DAS_stat{i_condition}.dv_mod(3)).^2);
id_temp = id_LT & (DAS_all{i_condition}.MovieNum==1) & (idx{i_condition} == 3);
DAS_track_id = DAS_all{i_condition}.TrackID(id_temp);
DAS_dist = dist_temp(id_temp);
Z1=Z_norm(id_temp,1);
Z2=Z_norm(id_temp,2);
[~, id_sort] = sort(DAS_dist);
id = DAS_track_id(id_sort);
id = id(1:num_track_select);
Z1 = Z1(id_sort);Z1=Z1(1:num_track_select);
Z2 = Z2(id_sort);Z2=Z2(1:num_track_select);
plotExampleTrack(data{i_condition}, Track_info,id,{gcf},pm,'subplot_id',[1,2,1:2]);
if pm.fig1_mosaic
subplot(1,2,2); 
end
xlim([0 30]); ylim([0 120]);
figure(figH{2}(2));

scatter(Z2,Z1,'filled','SizeData',size_exp,'MarkerFaceColor',col_scatter,'MarkerEdgeColor',col_scatter_edge);alpha(1);
%scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'+','MarkerEdgeColor',[0 1 0],'LineWidth',3);hold on;
%scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'x','MarkerEdgeColor',[0 1 0],'LineWidth',3);hold on;
scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'o','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;
scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'d','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;
%dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_var, DAS_all{i_condition}.DAS , idx{i_condition},DAS_stat{i_condition}.edges_d_dv, {label1,label2},'fig_exist',figH{2}(2),'contour_only',true);
%==========================================================================
figure(figH{2}(7));

dist_temp = sqrt((Z_norm(:,2)-DAS_stat{i_condition}.d_mod(1)).^2+(Z_norm(:,1)-DAS_stat{i_condition}.dv_mod(1)).^2);
id_temp = id_Imax & (DAS_all{i_condition}.MovieNum==1) & (idx{i_condition} == 1);
DAS_track_id = DAS_all{i_condition}.TrackID(id_temp);
DAS_dist = dist_temp(id_temp);
Z1=Z_norm(id_temp,1);
Z2=Z_norm(id_temp,2);
[~, id_sort] = sort(DAS_dist);
id = DAS_track_id(id_sort);
id = id(1:num_track_select);
Z1 = Z1(id_sort);Z1=Z1(1:num_track_select);
Z2 = Z2(id_sort);Z2=Z2(1:num_track_select);
plotExampleTrack(data{i_condition}, Track_info,id,{gcf},pm,'subplot_id',[1,2,1:2]);
if pm.fig1_mosaic
subplot(1,2,2); 
end
xlim([0 70]); ylim([0 120]);
figure(figH{2}(4));
dasClusterPlot_single(DAS_stat{i_condition}.z_d_dv, DAS_stat{i_condition}.edges_d_dv, {label2,label1},'fig_exist',figH{2}(4));
hold on
scatter(Z2,Z1,'filled','SizeData',size_exp,'MarkerFaceColor',col_scatter,'MarkerEdgeColor',col_scatter_edge);alpha(1);hold on;
%scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'+','MarkerEdgeColor',[0 1 0],'LineWidth',3);hold on;
%scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'x','MarkerEdgeColor',[0 1 0],'LineWidth',3);hold on;
scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'o','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;
scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'d','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;

figure(figH{2}(8));
dist_temp = sqrt((Z_norm(:,2)-DAS_stat{i_condition}.d_mod(3)).^2+(Z_norm(:,1)-DAS_stat{i_condition}.dv_mod(3)).^2);
id_temp = id_Imax & (DAS_all{i_condition}.MovieNum==1) & (idx{i_condition} == 3);
DAS_track_id = DAS_all{i_condition}.TrackID(id_temp);
DAS_dist = dist_temp(id_temp);
Z1=Z_norm(id_temp,1);
Z2=Z_norm(id_temp,2);
[~, id_sort] = sort(DAS_dist);
id = DAS_track_id(id_sort);
id = id(1:num_track_select);
Z1 = Z1(id_sort);Z1=Z1(1:num_track_select);
Z2 = Z2(id_sort);Z2=Z2(1:num_track_select);
plotExampleTrack(data{i_condition}, Track_info,id,{gcf},pm,'subplot_id',[1,2,1:2]);
if pm.fig1_mosaic
subplot(1,2,2); 
end
xlim([0 70]); ylim([0 120]);
figure(figH{2}(4));
scatter(Z2,Z1,'filled','SizeData',size_exp,'MarkerFaceColor',col_scatter,'MarkerEdgeColor',col_scatter_edge);alpha(1);
hold on
%dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_var, DAS_all{i_condition}.DAS , idx{i_condition},DAS_stat{i_condition}.edges_d_dv, {label1,label2},'fig_exist',figH{2}(4),'contour_only',true);


