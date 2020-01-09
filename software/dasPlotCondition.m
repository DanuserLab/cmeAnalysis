

function [figH] = dasPlotCondition(data, DAS_all, idx, con_name, DAS_stat, test_p,pm,dir_alt, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) iscell(x));
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('idx', @(x) iscell(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('DAS_stat', @(x) iscell(x));
ip.addRequired('test_p', @(x) isstruct(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addParameter('fig_name', 'fig_cond', @ischar);

ip.parse(data, DAS_all, idx, con_name, DAS_stat, test_p,pm,dir_alt, varargin{:});


DAS_all = ip.Results.DAS_all;
con_name = ip.Results.con_name;
idx = ip.Results.idx;
DAS_stat = ip.Results.DAS_stat;
fig_name = ip.Results.fig_name;
pm = ip.Results.pm;
dir_alt = ip.Results.dir_alt;
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
num_condition = max(size(DAS_all));
num_clus =3;
%==========================================================================
label1 = pm.label1;
label2 = pm.label2;
label3 = pm.label3;
%==========================================================================
color_clus = pm.color_clus;
%==========================================================================
figH = cell(3,1);
%==========================================================================
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.addGroup(fig_name);
desktop.setGroupDocked(fig_name, 0);
myDim   = java.awt.Dimension(num_condition, 3);   % 4 columns, 2 rows
% 1: Maximized, 2: Tiled, 3: Floating
n_fig = 3 * num_condition;
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
for i_condition = 1:num_condition
    figure(figH{1}(i_condition)); hold on
    title(con_name{i_condition});
dasClusterPlot_single(DAS_stat{i_condition}.z_d_dv, DAS_stat{i_condition}.edges_d_dv, {label2,label1},'fig_exist',figH{1}(i_condition));
if pm.plot_mode
scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'o','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;
scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'x','MarkerEdgeColor',[0 1 0],'LineWidth',1);hold on;
end
figure(figH{1}(i_condition+num_condition)); hold on
dasClusterPlot_single_contour(DAS_all{i_condition}.DAS_var, DAS_all{i_condition}.DAS , idx{i_condition},DAS_stat{i_condition}.edges_d_dv, {label2,label1},pm,'fig_exist',figH{1}(i_condition+num_condition));
if pm.plot_mode
scatter(DAS_stat{i_condition}.d_mod(1), DAS_stat{i_condition}.dv_mod(1),100,'o','MarkerEdgeColor',[1 1 1],'LineWidth',1);hold on;
scatter(DAS_stat{i_condition}.d_mod(3), DAS_stat{i_condition}.dv_mod(3),100,'x','MarkerEdgeColor',[0 1 0],'LineWidth',1);hold on;
end
if i_condition > 1
figure(figH{1}(i_condition+num_condition*2));
sum_tem1 = sum(sum(DAS_stat{i_condition}.z_d_dv));
sum_tem2 = sum(sum(DAS_stat{1}.z_d_dv));
Delta_d1d2=(DAS_stat{1}.edges_d_dv{2}(end)-DAS_stat{1}.edges_d_dv{2}(end-1))...
          *(DAS_stat{1}.edges_d_dv{1}(end)-DAS_stat{1}.edges_d_dv{1}(end-1));
dasClusterPlot_single(DAS_stat{i_condition}.z_d_dv/sum_tem1/Delta_d1d2-DAS_stat{1}.z_d_dv/sum_tem2/Delta_d1d2,...
          DAS_stat{i_condition}.edges_d_dv, {label2,label1},...
          'fig_exist',figH{1}(i_condition+num_condition*2),'is_diff',true...
          );
end
end
%==========================================================================
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.addGroup([fig_name '_2']);
desktop.setGroupDocked([fig_name '_2'], 0);
myDim   = java.awt.Dimension(4, 3);   % 4 columns, 2 rows
% 1: Maximized, 2: Tiled, 3: Floating
n_fig = 4 * 3;
figH{2}    = gobjects(n_fig, 1);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
for iFig = 1:n_fig
    if (iFig > 6) && (iFig<11) 
        if  (pm.pdf_conf == false)
        figH{2}(iFig) = figure('WindowStyle', 'docked', ...
         'Name', sprintf('Figure %d', iFig), 'NumberTitle', 'off');
        set(get(handle(figH{2}(iFig)), 'javaframe'), 'GroupName', [fig_name '_2']);
        end
    else
    figH{2}(iFig) = figure('WindowStyle', 'docked', ...
      'Name', sprintf('Figure %d', iFig), 'NumberTitle', 'off');
    set(get(handle(figH{2}(iFig)), 'javaframe'), 'GroupName', [fig_name '_2']);
    end
end
if pm.pdf_conf == true
   fig_exist = dasPdfConf(data, dir_alt, pm, con_name, DAS_stat, DAS_all, idx,'LT_not_Imax',true);
   figH{2}(7) = fig_exist{1}; figH{2}(8) = fig_exist{2};
   set(figH{2}(7),'WindowStyle','docked');
   set(get(handle(figH{2}(7)), 'javaframe'), 'GroupName', [fig_name '_2']);
   %---------------------------
   set(figH{2}(8),'WindowStyle','docked');
   set(get(handle(figH{2}(8)), 'javaframe'), 'GroupName', [fig_name '_2']);
   %---------------------------
   fig_exist = dasPdfConf(data, dir_alt, pm, con_name, DAS_stat, DAS_all, idx,'LT_not_Imax',false);
   figH{2}(9) = fig_exist{1}; figH{2}(10) = fig_exist{2};
   set(figH{2}(9),'WindowStyle','docked');
   set(get(handle(figH{2}(9)), 'javaframe'), 'GroupName', [fig_name '_2']);
   %---------------------------
   set(figH{2}(10),'WindowStyle','docked');
   set(get(handle(figH{2}(10)), 'javaframe'), 'GroupName', [fig_name '_2']);
   %---------------------------
end
warning(bakWarn);
desktop.setDocumentArrangement([fig_name '_2'], 2, myDim)
%==========================================================================
dasPlotBar(DAS_all, idx, DAS_stat,test_p,con_name, pm, 'fig_exist',{figH{2}(1),figH{2}(2),figH{2}(3),figH{2}(4),figH{2}(5)});
%==========================================================================
figure(figH{2}(6));
x_temp = [];
for i_condition = 1:num_condition
    if pm.plot_visitor == true
        x_temp = cat(1,x_temp,DAS_stat{i_condition}.prop_c');
    else
        x_temp = cat(1,x_temp,[DAS_stat{i_condition}.prop_c(1)/(DAS_stat{i_condition}.prop_c(1)+DAS_stat{i_condition}.prop_c(3))...
                               DAS_stat{i_condition}.prop_c(3)/(DAS_stat{i_condition}.prop_c(1)+DAS_stat{i_condition}.prop_c(3))]);
    end
end
if num_condition > 1
    c = categorical(con_name);
    c = reordercats(c,con_name);
h = bar(c,x_temp,'stacked');
if pm.plot_visitor == true
set(h,{'FaceColor'},{color_clus{1};color_clus{2};color_clus{3}});
else
    set(h,{'FaceColor'},{color_clus{1};color_clus{3}});
end
else
    if pm.plot_visitor == true
       c = categorical({'CCPs','FPs','PCs'});
    else
        c = categorical({'CCP','PCs'});
    end
   h = bar(c, x_temp);
end
xtickangle(45)
%==========================================================================
if pm.pdf_conf == false
%==========================================================================
figure(figH{2}(7));
title('Lifetime Distribution of CCPs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.LT1,DAS_stat{i_condition}.pdf_LT1,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
ylabel('Prob. Density');xlabel('$\tau$ (s)','interpreter','latex')
xlim([0 150])
end
%==========================================================================
figure(figH{2}(8));
title('I_{max} Distribution of CCPs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.Imax1,DAS_stat{i_condition}.pdf_Imax1,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
end
%legend(con_name)
ylabel('Prob. Density');xlabel('$I_{max}$ (a.u.)','interpreter','latex')
xlim([0 250])
%==========================================================================
figure(figH{2}(9));
title('Lifetime Distribution of PCs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.LT3,DAS_stat{i_condition}.pdf_LT3,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
end
%legend(con_name)
ylabel('Prob. Density');xlabel('$\tau$ (s)','interpreter','latex')
xlim([0 50])
%==========================================================================
figure(figH{2}(10));
title('I_{max} Distribution of PCs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.Imax3,DAS_stat{i_condition}.pdf_Imax3,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
end
%legend(con_name)
ylabel('Prob. Density');xlabel('$I_{max}$ (a.u.)','interpreter','latex')
xlim([0 250])
%==========================================================================
end
%==========================================================================
figure(figH{2}(11));
title('Lifetime Distribution of FPs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.LT2,DAS_stat{i_condition}.pdf_LT2,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
end
%legend(con_name)
ylabel('Prob. Density');xlabel('$\tau$ (s)','interpreter','latex')
xlim([0 150])
%==========================================================================
figure(figH{2}(12));
title('I_{max} Distribution of FPs')
hold on
for i_condition = 1:num_condition
line(DAS_stat{i_condition}.Imax2,DAS_stat{i_condition}.pdf_Imax2,'Color',pm.col_cond(i_condition,:),'Linewidth',1)
end
%legend(con_name)
ylabel('Prob. Density');xlabel('$I_{max}$ (a.u.)','interpreter','latex')
xlim([0 250])
%==========================================================================
% %write signifiance in figures
% %---------------------------
%    figure(figH{2}(7))
%    xtem = get(gca,'XLim'); ytem = get(gca,'YLim'); i_c = 1;
%    if pm.show_med_or_mean; txt_tem = {'median \tau '}; else;  txt_tem = {'mean \tau '}; end
%     for i_condition = 2:num_condition
%      if pm.show_med_or_mean; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_LT_med{i_c});
%      else; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_LT_mean{i_c});
%      end
%     end
%    text(xtem(1)+(xtem(2)-xtem(1))*0.5,ytem(1)+(ytem(2)-ytem(1))*0.8,txt_tem,'FontSize',8)
%    %---------------------------
%    figure(figH{2}(8))
%    xtem = get(gca,'XLim'); ytem = get(gca,'YLim'); i_c = 3;
%    if pm.show_med_or_mean; txt_tem = {'median \tau '}; else;  txt_tem = {'mean \tau '}; end
%     for i_condition = 2:num_condition
%      if pm.show_med_or_mean; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_LT_med{i_c});
%      else; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_LT_mean{i_c});
%      end
%     end
%    text(xtem(1)+(xtem(2)-xtem(1))*0.5,ytem(1)+(ytem(2)-ytem(1))*0.8,txt_tem,'FontSize',8)
%    %---------------------------
%    figure(figH{2}(9))
%    xtem = get(gca,'XLim'); ytem = get(gca,'YLim'); i_c = 1;
%    if pm.show_med_or_mean; txt_tem = {'median{\it I_{max}} '}; else;  txt_tem = {'mean{\it I_{max}} '}; end
%     for i_condition = 2:num_condition
%      if pm.show_med_or_mean; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_Imax_med{i_c});
%      else; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_Imax_mean{i_c});
%      end
%     end
%    text(xtem(1)+(xtem(2)-xtem(1))*0.5,ytem(1)+(ytem(2)-ytem(1))*0.8,txt_tem,'FontSize',8)
%    %---------------------------
%    figure(figH{2}(10))
%    xtem = get(gca,'XLim'); ytem = get(gca,'YLim'); i_c = 3;
%    if pm.show_med_or_mean; txt_tem = {'median{\it I_{max}} '}; else;  txt_tem = {'mean{\it I_{max}} '}; end
%     for i_condition = 2:num_condition
%      if pm.show_med_or_mean; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_Imax_med{i_c});
%      else; txt_tem = cat(1,txt_tem,DAS_stat{i_condition}.ps_Imax_mean{i_c});
%      end
%     end
%    text(xtem(1)+(xtem(2)-xtem(1))*0.5,ytem(1)+(ytem(2)-ytem(1))*0.8,txt_tem,'FontSize',8)
%    %---------------------------
%==========================================================================
if pm.plot_cohort
[fig_exist] = dasCohorts(data, DAS_all, idx,dir_alt, pm);
%==========================================================================
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.addGroup([fig_name '_3']);
desktop.setGroupDocked([fig_name '_3'], 0);
if pm.plot_EpiTIRF == true
myDim   = java.awt.Dimension(size(fig_exist,1)+2, size(fig_exist,2));   % 4 columns, 2 rows
n_fig = (size(fig_exist,2))* (size(fig_exist,1)+2);
else
    myDim   = java.awt.Dimension(size(fig_exist,1), size(fig_exist,2));
    n_fig = size(fig_exist,2)* (size(fig_exist,1));
end

% 1: Maximized, 2: Tiled, 3: Floating

figH{3}    = gobjects(n_fig, 1);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
for iFig = 1:size(fig_exist,2)* (size(fig_exist,1))
   figH{3}(iFig) = fig_exist{iFig};
   set(figH{3}(iFig),'WindowStyle','docked');
   set(get(handle(figH{3}(iFig)), 'javaframe'), 'GroupName', [fig_name '_3']);
end
warning(bakWarn);
%==========================================================================
if pm.plot_EpiTIRF == true   
    fig_exist = dasEpiTIRF(data, dir_alt, pm, con_name, DAS_stat,DAS_all,idx);
    figH{3}(end-1) = fig_exist{1};
    figH{3}(end) = fig_exist{2};
    set(figH{3}(end-1),'WindowStyle','docked');
    set(figH{3}(end),'WindowStyle','docked');
    set(get(handle(figH{3}(end-1)), 'javaframe'), 'GroupName', [fig_name '_3']);
    set(get(handle(figH{3}(end)), 'javaframe'), 'GroupName', [fig_name '_3']);
end

desktop.setDocumentArrangement([fig_name '_3'], 2, myDim)
end
