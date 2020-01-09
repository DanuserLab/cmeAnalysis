

function [figH] = dasFigure1_3D(data, DAS_all,Track_info, DAS_stat, idx, Cx, con_name,dir_alt,pm,figH,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) iscell(x));
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addRequired('DAS_stat', @(x) iscell(x));
ip.addRequired('idx', @(x) iscell(x));
ip.addRequired('Cx', @(x) iscell(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('figH', @(x) iscell(x));
ip.addParameter('fix_cluster_num', 2, @isnumeric);
ip.addParameter('Z_rep', [-0.2, -4; -0.8, -2], @isnumeric);
ip.addParameter('i_condition', 1, @isnumeric);
ip.addParameter('fig_name', 'fig1_3D', @ischar);
ip.parse(data, DAS_all, Track_info,DAS_stat, idx, Cx, con_name,dir_alt, pm,figH,varargin{:});


DAS_all = ip.Results.DAS_all;
con_name = ip.Results.con_name;
num_condition = max(size(DAS_all));
idx = ip.Results.idx;
Z_rep = ip.Results.Z_rep;
data = ip.Results.data;
dir_alt = ip.Results.dir_alt;
Track_info = ip.Results.Track_info;
DAS_stat = ip.Results.DAS_stat;
Cx = ip.Results.Cx;
i_condition = ip.Results.i_condition;
fig_name = ip.Results.fig_name;
pm=ip.Results.pm;
figH = ip.Results.figH;
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

%==========================================================================
    X = DAS_all{i_condition}.DAS_var;
    Y = DAS_all{i_condition}.DAS;
    Z = DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3;
    Z_org=[X,Y,Z];
    Z_norm=[normalize(X,'zscore'),normalize(Y,'zscore'),normalize(Z,'zscore')];
    %-------------------------------------------------------------------------
    axis_lim = [min(DAS_stat{i_condition}.edges_d_dv{1}) max(DAS_stat{i_condition}.edges_d_dv{1});...
                min(DAS_stat{i_condition}.edges_d_dv{2}) max(DAS_stat{i_condition}.edges_d_dv{2});...
                min(DAS_stat{i_condition}.edges_d_d32{1}) max(DAS_stat{i_condition}.edges_d_d32{1})];
    r = 0.2*(axis_lim(:,2)-axis_lim(:,1))/sqrt(sum((axis_lim(:,2)-axis_lim(:,1)).^2));  
    %-------------------------------------------------------------------------
%     Cx_org = Cx;
%     Cx_org{i_condition}(:,1) = Cx{i_condition}(:,1).*std(X)+mean(X);
%     Cx_org{i_condition}(:,2) = Cx{i_condition}(:,2).*std(Y)+mean(Y);
%     Cx_org{i_condition}(:,3) = Cx{i_condition}(:,3).*std(Z)+mean(Z);

    %-------------------------------------------------------------------------
    %ZC_dist = [sum((Z_org - Cx_org{i_condition}(1,:)).^2,2), sum((Z_org - Cx_org{i_condition}(2,:)).^2,2), sum((Z_org - Cx_org{i_condition}(3,:)).^2,2)];
    %d_temp = 0.02;
    %-------------------------------------------------------------------------
%     id_temp = (abs(ZC_dist(:,1)-ZC_dist(:,2)) < d_temp) & ((idx{i_condition} == 1) | (idx{i_condition} == 2)); 
%     X_temp = [ones(size(Z_org(id_temp,1))) Z_org(id_temp,1) Z_org(id_temp,2) Z_org(id_temp,1).*Z_org(id_temp,2)];
%     b_temp = regress(Z_org(id_temp,3),X_temp) ;
%     x1fit = min(Z_org(id_temp,1)):0.1:max(Z_org(id_temp,1));
%     x2fit = min(Z_org(id_temp,2)):0.1:max(Z_org(id_temp,2));
%     [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
%     YFIT = b_temp(1) + b_temp(2)*X1FIT + b_temp(3)*X2FIT + b_temp(4)*X1FIT.*X2FIT;
%     h_temp = mesh(X1FIT,X2FIT,YFIT);
%     set(h_temp,'EdgeColor','none','FaceColor', [0 1 0]);   
%     hold on
%     %-------------------------------------------------------------------------
%     id_temp = (abs(ZC_dist(:,1)-ZC_dist(:,3)) < d_temp) & ((idx{i_condition} == 1) | (idx{i_condition} == 3)); 
%     X_temp = [ones(size(Z_org(id_temp,1))) Z_org(id_temp,1) Z_org(id_temp,2) Z_org(id_temp,1).*Z_org(id_temp,2)];
%     b_temp = regress(Z_org(id_temp,3),X_temp) ;
%     x1fit = min(Z_org(id_temp,1)):0.1:max(Z_org(id_temp,1));
%     x2fit = min(Z_org(id_temp,2)):0.1:max(Z_org(id_temp,2));
%     [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
%     YFIT = b_temp(1) + b_temp(2)*X1FIT + b_temp(3)*X2FIT + b_temp(4)*X1FIT.*X2FIT;
%     h_temp = mesh(X1FIT,X2FIT,YFIT);
%     set(h_temp,'EdgeColor','none','FaceColor', [0 0 1]);   
%     hold on
    %-------------------------------------------------------------------------
%     id_temp_all = ((abs(ZC_dist(:,1)-ZC_dist(:,3) )< d_temp )) & ((idx{i_condition} == 1) | (idx{i_condition} == 3)); 
%     f_temp = fit([Z_org(id_temp_all,1), Z_org(id_temp_all,2)],Z_org(id_temp_all,3),'linearinterp');
%     h_temp=plot(f_temp);
%     set(h_temp,'EdgeColor','none','FaceAlpha',1, 'FaceColor',color_clus{3});
%     xlim(axis_lim(1,:));
%     ylim(axis_lim(2,:));   
%     zlim(axis_lim(3,:));
%     hold on
%     id_temp_all = ((abs(ZC_dist(:,1)-ZC_dist(:,2) )< d_temp )) & ((idx{i_condition} == 1) | (idx{i_condition} == 2)); 
%     f_temp = fit([Z_org(id_temp_all,1), Z_org(id_temp_all,2)],Z_org(id_temp_all,3),'linearinterp');   
%     h_temp=plot(f_temp);
%     set(h_temp,'EdgeColor','none','FaceAlpha',1, 'FaceColor',color_clus{2});
%     hold on
%-------------------------------------------------------------------------

%==========================================================================
    %-------------------------------------------------------------------------
    figure(figH{1}(1));
    %-------------------------------------------------------------------------
    alpha_tem = 0.8;
colormap('bone');
[X,Y] = meshgrid(DAS_stat{i_condition}.edges_d_dv{1},DAS_stat{i_condition}.edges_d_dv{2});
z=transpose(DAS_stat{i_condition}.z_d_dv/max(max(DAS_stat{i_condition}.z_d_dv)));
hSurface = surf(X, Y,z+axis_lim(3,1),z);
set(hSurface,'FaceLighting','gouraud','EdgeColor','none','FaceColor', 'interp','FaceAlpha',alpha_tem);
%mesh(X, Y,z+axis_lim(3,1),z,'FaceLighting','gouraud','LineWidth',1.5)

hold on

[X,Y] = meshgrid(DAS_stat{i_condition}.edges_d_d32{1},DAS_stat{i_condition}.edges_d_d32{2});
z=transpose(DAS_stat{i_condition}.z_d_d32/max(max(DAS_stat{i_condition}.z_d_d32)));
hSurface = surf(X, Y,z,z);
direction = [0 1 0];
rotate(hSurface,direction,-90)
trans = max(max(get(hSurface,'XData')));
XData_temp = get(hSurface,'XData');
set(hSurface,'XData',XData_temp-trans+axis_lim(1,2),'FaceLighting','gouraud','EdgeColor','none','FaceColor', 'interp','FaceAlpha',alpha_tem);
hold on

[X,Y] = meshgrid(DAS_stat{i_condition}.edges_d32_dv{1},DAS_stat{i_condition}.edges_d32_dv{2});
z=transpose(DAS_stat{i_condition}.z_d32_dv/max(max(DAS_stat{i_condition}.z_d32_dv)));
hSurface = surf(X, Y,z,z);
direction = [1 0 0];
rotate(hSurface,direction,90)
trans = max(max(get(hSurface,'YData')));
YData_temp = get(hSurface,'YData');
set(hSurface,'YData',YData_temp-trans+axis_lim(2,2),'FaceLighting','gouraud','EdgeColor','none','FaceColor', 'interp','FaceAlpha',alpha_tem);
hold on

colorbar 
caxis([0 1])
    xlim(axis_lim(1,:));
    ylim(axis_lim(2,:));   
    zlim(axis_lim(3,:));
c = colorbar;
c.Ticks = [];
grid off;
%========================================================================== 
 
%========================================================================== 
    %-------------------------------------------------------------------------
    figure(figH{1}(1));
    %-------------------------------------------------------------------------
for i_c = 1:num_clus
    num_track_select = floor(200*DAS_stat{i_condition}.prop_c(i_c)+0.5);
dist_temp = sqrt((Z_org(:,2)-DAS_stat{i_condition}.d_mod(i_c)).^2+(Z_org(:,1)-DAS_stat{i_condition}.dv_mod(i_c)).^2+(Z_org(:,3)-DAS_stat{i_condition}.d32_mod(i_c)).^2);
id_temp = (idx{i_condition} == i_c);
DAS_dist = dist_temp(id_temp);
Z1=Z_org(id_temp,1);
Z2=Z_org(id_temp,2);
Z3=Z_org(id_temp,3);
%[~, id_sort] = sort(DAS_dist);
id_sort = randsample(max(size(DAS_dist)),num_track_select);
Z1 = Z1(id_sort);Z1=Z1(1:num_track_select);
Z2 = Z2(id_sort);Z2=Z2(1:num_track_select);
Z3 = Z3(id_sort);Z3=Z3(1:num_track_select);

for i_plot = 1:num_track_select
[x,y,z] = sphere(50);
s1 = surf(x*r(1)+Z1(i_plot),y*r(2)+Z2(i_plot),z*r(3)+Z3(i_plot));
hold on

%axis equal

set(s1,'facecolor',color_clus{i_c},'FaceLighting','gouraud','EdgeColor','none','facealpha',1); 
%xlim([-4 4]);ylim([-4 4]);zlim([-4 4]);
end
hold on
end

xlabel(label1,'interpreter','latex');
ylabel(label2,'interpreter','latex');
zlabel(label3,'interpreter','latex');

%light('Position',[-0.5 0.2 0.9],'Style','infinite');
light('Position',[-0.2 -0.2 1],'Style','infinite');
%camlight
view([-30 30]);










