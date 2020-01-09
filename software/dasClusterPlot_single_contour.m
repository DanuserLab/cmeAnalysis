

function [] = dasClusterPlot_single_contour(X, Y , idx,edges, var_name,pm,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('X', @(x) isnumeric(x));
ip.addRequired('Y', @(x) isnumeric(x));
ip.addRequired('idx', @(x) isnumeric(x));
ip.addRequired('edges', @(x) iscell(x));
ip.addRequired('var_name', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('fix_cluster_num', 2, @isnumeric);
ip.addParameter('fig_exist', [], @ishandle);
ip.addParameter('contour_only', false, @islogical);
ip.addParameter('exp_only', true, @islogical);
ip.parse(X, Y, idx,edges, var_name, pm,varargin{:});


X = ip.Results.X;
Y = ip.Results.Y;
var_name = ip.Results.var_name;
pm = ip.Results.pm;
idx = ip.Results.idx;
edges = ip.Results.edges;
fig_exist = ip.Results.fig_exist;
contour_only = ip.Results.contour_only;
exp_only = ip.Results.exp_only;
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
num_clus = 3;

color_clus = pm.color_clus;
%==========================================================================
if isempty(fig_exist)
    fig_cluster = figure;
    %set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 10, 20]);
else
    fig_cluster = fig_exist;
end

figure(fig_cluster);
    Z_norm=[X,Y];
  %---------------------------------------------------- 
   if contour_only == false
      for i_c = 1: num_clus
          id_temp = (idx==i_c);
          x_temp = Z_norm(id_temp,1:2);
        x_temp = randSampling(x_temp,0.1);     
        scatter(x_temp(:,2),x_temp(:,1),'filled','SizeData',5,'MarkerFaceColor',color_clus{i_c});alpha(.1);
        hold on;
      end
   end
   if exp_only == false
      z = hist3(Z_norm,'Edges',edges);
      [~,ct] = contour(edges{2},edges{1},z/max(max(z)),v_cont,'--','LineWidth', 0.5);
      colorbar; caxis([0 1])
   end
      xlabel(var_name{1},'interpreter','latex'); ylabel(var_name{2},'interpreter','latex'); 
   xlim([min(edges{2}) max(edges{2})]);ylim([min(edges{1}) max(edges{1})]);  

   %---------------------------------------------------- 