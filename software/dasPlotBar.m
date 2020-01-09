

function [test_target,fig_bar] = dasPlotBar(DAS_all, idx, DAS_stat, test_p, con_name, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('idx', @(x) iscell(x));
ip.addRequired('DAS_stat', @(x) iscell(x));
ip.addRequired('test_p', @(x) isstruct(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('fix_cluster_num', 2, @isnumeric);
ip.addParameter('fig_exist', [], @iscell);
ip.parse(DAS_all, idx, DAS_stat, test_p, con_name, pm, varargin{:});


DAS_all = ip.Results.DAS_all;
con_name = ip.Results.con_name;
idx = ip.Results.idx;
num_condition = size(DAS_all,1);
fig_exist = ip.Results.fig_exist;
pm = ip.Results.pm;
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
n_bar = pm.n_bar;

fig_bar = cell(1,1);
fig_bar{1}    = gobjects(n_bar, 1);
if isempty(fig_exist)
    for i_bar = 1:n_bar
        fig_bar{1}(i_bar) = figure;
    end
else
    for i_bar = 1:n_bar
    fig_bar{1}(i_bar) = fig_exist{i_bar};
    end
end


test_target = test_p.test_target;

p = zeros(num_condition-1,1);
pair = cell(num_condition-1,1);
for i_bar = 1:n_bar
    figure(fig_bar{1}(i_bar));
        if i_bar == 1
            ylabel('CCP (min^{-1}{\mu}m^{-2})')
            ylim([0 0.3])
        elseif i_bar == 2
            ylabel('OT (min^{-1}{\mu} m^{-2})')
            ylim([0 0.2])
        elseif i_bar == 3
            ylabel('AP (min^{-1}{\mu}m^{-2})')
            ylim([0 0.4])
        elseif i_bar == 4
            ylabel('CS initiation (min^{-1}{\mu}m^{-2})')
            ylim([0 1.1])
        elseif i_bar == 5
            ylabel('CCP%')
            ylim([0 60])
        end
        %----------------------------------
        for i_condition = 2:num_condition
            if (pm.scheme_p == 1) || (i_condition == 2)
                p(i_condition-1) = ranksum(test_target{i_bar}{1},test_target{i_bar}{i_condition});
            elseif (pm.scheme_p == 2) && (i_condition > 2)
                p(i_condition-1) = ranksum(test_target{i_bar}{2},test_target{i_bar}{i_condition});
            end
        end
        if i_bar < n_bar+1
        %----------------------------------
        h=notBoxPlot(test_target{i_bar}{1},1);
	    set([h.data],'MarkerSize',1)
        hold on

        for i_condition = 2:num_condition
        h=notBoxPlot(test_target{i_bar}{i_condition},i_condition);
	    set([h.data],'MarkerSize',1) 
        hold on
        end
        
        %----------------------------------
        for i_condition = 2:num_condition
        if (p(i_condition-1)<=0.05) &&  (p(i_condition-1)>1E-2)
          p(i_condition-1) = 0.05;
        elseif (p(i_condition-1)<=1E-2) &&  (p(i_condition-1)>1E-3)
          p(i_condition-1) = 1E-2;
        elseif (p(i_condition-1)<=1E-3) 
          p(i_condition-1) = 1E-3;
        else
          p(i_condition-1) = nan;
        end
        if (pm.scheme_p == 1) || (i_condition == 2)
            pair{i_condition-1} = [1,i_condition];
        elseif (pm.scheme_p == 2) && (i_condition > 2)
            pair{i_condition-1} = [2,i_condition];
        end
        end
        %----------------------------------
	    sigstar(pair,p,1);
	    box on
        xticks((1:num_condition))
        set(gca,'XTickLabel',con_name)
        xtickangle(45)
        %----------------------------------
        end
     %set(gca,'fontsize',17, 'FontWeight', 'bold')
     %set(gca,'Linewidth',2) 
end