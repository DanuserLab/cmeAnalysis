

function [figH] = dasSavePlot(figH, dir_alt, i_fig, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('figH', @(x) iscell(x));
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('i_fig', @(x) isnumeric(x));
ip.addParameter('PaperPosition', [0 0 5 5], @isnumeric);

ip.parse(figH, dir_alt, i_fig, varargin{:});


figH = ip.Results.figH;
dir_alt = ip.Results.dir_alt;
i_fig = ip.Results.i_fig;
PaperPosition = ip.Results.PaperPosition;
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
n_g = max(size(figH));
for i_g = 1: n_g
    dir_temp = [dir_alt filesep 'fig' num2str(i_fig) filesep 'group' num2str(i_g)];
    mkdir(dir_temp);
    n_f = max(size(figH{i_g}));
    for i_f = 1: n_f
        figure(figH{i_g}(i_f));
      set(figH{i_g}(i_f),'WindowStyle','normal');
%       gca_all = get(gcf,'Children');
%       for i_gca = 1: max(size(gca_all))
%           set(gca(i_gca),'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
%       end
    set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
      
      lgd=findobj(gcf,'type','legend');
      if ~isempty(lgd)
      legend('boxoff');
      legend('Location','northeast');
      end
      
      %set(gca,'fontsize',8, 'FontWeight', 'bold');
      %set(figH{i_g}(i_f), 'Units', 'Centimeters', 'OuterPosition', [0, 0, 5, 5]);   
      %saveas(figH{i_g}(i_f),[dir_temp filesep 'fig_' num2str(i_fig) '_g_' num2str(i_g) '_no_' num2str(i_f)],'png');	  
	  %figH{i_g}(i_f).PaperPositionMode = 'auto';
      figH{i_g}(i_f).PaperUnits = 'centimeters';
      figH{i_g}(i_f).PaperPosition = PaperPosition;
	  print(figH{i_g}(i_f),[dir_temp filesep 'fig_' num2str(i_fig) '_g_' num2str(i_g) '_no_' num2str(i_f)],'-dpng','-r900')
      set(figH{i_g}(i_f),'WindowStyle','docked');
      savefig(figH{i_g}(i_f),[dir_temp filesep 'fig_' num2str(i_fig) '_g_' num2str(i_g) '_no_' num2str(i_f)])
    end

end
