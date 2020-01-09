

function [] = dasClusterPlot_single(z , edges, var_name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('z', @(x) isnumeric(x));
ip.addRequired('edges', @(x) iscell(x));
ip.addRequired('var_name', @(x) iscell(x));
ip.addParameter('fix_cluster_num', 2, @isnumeric);
ip.addParameter('fig_exist', [], @ishandle);
ip.addParameter('map_col', 'jet', @ischar);
ip.addParameter('is_diff', false, @islogical);
ip.addParameter('alt_norm', 1, @isnumeric);
ip.parse(z, edges, var_name, varargin{:});


z = ip.Results.z;
var_name = ip.Results.var_name;

edges = ip.Results.edges;
fig_exist = ip.Results.fig_exist;
map_col = ip.Results.map_col;
is_diff = ip.Results.is_diff;
alt_norm = ip.Results.alt_norm;
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
v_cont = (0:0.1:1);
if is_diff
    v_cont = (-0.4:0.05:0.4);
end
%==========================================================================
if isempty(fig_exist)
    fig_cluster = figure;
    %set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 10, 20]);
else
    fig_cluster = fig_exist;
end

figure(fig_cluster)
  %----------------------------------------------------  
   %colormap('hot');
   if is_diff
       colormap('hot');
   else
       colormap(map_col);
   end
   if is_diff
       z_tem = z/alt_norm;
   else
       z_tem = z/max(max(z))/alt_norm;
   end
   if is_diff
       z_tem(1) = min(v_cont);
       z_tem(end) = max(v_cont);
       contourf(edges{2},edges{1},z_tem,v_cont,':','LineWidth', 0.1);
   else
       contourf(edges{2},edges{1},z_tem,v_cont,'-','LineWidth', 1);
   end
   if is_diff
   xlim([min(edges{2}(2:end-1)) max(edges{2}(2:end-1))]);ylim([min(edges{1}(2:end-1)) max(edges{1}(2:end-1))]);  
   set(gca,'XTickLabel',[]);set(gca,'xtick',[]);
   set(gca,'YTickLabel',[]);set(gca,'ytick',[]);
   ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
   else
   colorbar; xlabel(var_name{1},'interpreter','latex'); ylabel(var_name{2},'interpreter','latex'); 
   xlim([min(edges{2}) max(edges{2})]);ylim([min(edges{1}) max(edges{1})]);
   end
  %---------------------------------------------------- 