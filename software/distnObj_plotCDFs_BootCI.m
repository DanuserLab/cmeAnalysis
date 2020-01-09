function figout = distnObj_plotCDFs_BootCI(DO_ctrl, DO_trt, ctlName, trtName, varName, varargin)
% distnObj_plotPDFs_BootCI VISUALIZE multiple distributional data objects 
% through cumulative distribution functions (CDFs) with bootstrapped confidence
% intervals.
% 
% Usage:
%   distnObj_plotCDFs_BootCI(DO_ctrl, DO_trt, 'condition1', 'condition2', ...
%       'TelomereLength', 'PositiveSupport', true)
%
% Input: 
%   - DO_ctrl is a struct with distributional data objects which can be obtained by
%   converting a pooled table using pooledTableToDistnObjects().
%
% Options:
%   - numBoots: default is 400.
%   - save: Default is true. If false, output is not saved.
%
% Updates:
% Jungsik Noh, 2018/02/13. Add 'save' option.
% Jungsik Noh, 2018/02/04.
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

ip = inputParser;
ip.addParameter('save', true);
ip.addParameter('outputDir', fullfile(pwd, 'distnObj_CDFs_BootCI'));
ip.addParameter('numBoots', 400)
ip.parse(varargin{:});
p = ip.Results;
 
%%

N1 = numel(DO_ctrl);
N2 = numel(DO_trt);

figout = cell(10,1);
  

%% 1. emprirical cdf

figout{1} = figure;
%[f, x] = ecdf(DASdistns{1}, 'function', 'cumulative hazard'); plot(x, f)
%ecdf(DASdistns{1})
xx =  DO_ctrl(1).uniqueX;
cdfval = DO_ctrl(1).cdfval;
stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])
grid on
ax = gca; ax.YLim=[0, 1];

hold on
for k = 2:numel(DO_ctrl)
    %[f, x] = ecdf(DASdistns{k}); plot(x, f)
    %ecdf(DASdistns{k})
    xx =  DO_ctrl(k).uniqueX;
    cdfval = DO_ctrl(k).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])

end

xlabel(varName)
ylabel('Cumulative frequency')
title([varName, '-', ctlName])
 

%
[xx, cdfval] = poolingDO(DO_ctrl);
stairs([xx(1); xx], [0; cdfval], 'Color', 'r', 'LineWidth', 2)

%% ecdf trt


figout{2} = figure;
%[f, x] = ecdf(DASdistns{1}, 'function', 'cumulative hazard'); plot(x, f)
%ecdf(DASdistns{1})
xx =  DO_trt(1).uniqueX;
cdfval = DO_trt(1).cdfval;
stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])
grid on
ax = gca; ax.YLim=[0, 1];

hold on
for k = 2:numel(DO_trt)
    %[f, x] = ecdf(DASdistns{k}); plot(x, f)
    %ecdf(DASdistns{k})
    xx =  DO_trt(k).uniqueX;
    cdfval = DO_trt(k).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])

end

xlabel(varName)
ylabel('Cumulative frequency')
title([varName, '-', trtName])
 
%
[xx, cdfval] = poolingDO(DO_trt);
stairs([xx(1); xx], [0; cdfval], 'Color', 'r', 'LineWidth', 2)


%% ecdfs together

colOrd = get(gca, 'ColorOrder');

figout{3} = figure;
    xx = DO_ctrl(1).uniqueX;
    cdfval = DO_ctrl(1).cdfval;
    p1 = stairs([xx(1); xx], [0; cdfval], 'Color', colOrd(1, :));

hold on

for i=2:N1
    xx = DO_ctrl(i).uniqueX;
    cdfval = DO_ctrl(i).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', colOrd(1, :))
    ylim([0, 1])
end

    xx = DO_trt(1).uniqueX;
    cdfval = DO_trt(1).cdfval;
    p2 = stairs([xx(1); xx], [0; cdfval], 'Color', colOrd(2, :));
    
for i=2:N2
    xx = DO_trt(i).uniqueX;
    cdfval = DO_trt(i).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', colOrd(2, :))
    ylim([0, 1])
end

legend([p1 p2], ctlName, trtName, 'Location', 'northeast')


%%
%%  2. Bootstrap Confidence Interval of CDFs
%%

DO = DO_ctrl;

% DO_ctrl
N = numel(DO);
[xall, fullcdf] = poolingDO(DO);

minx = min(xall); maxx = max(xall);
xgrid = linspace(minx, maxx, 1001);              % num of grid is set to 1001.
numel(xgrid);

%tic
cdfgridFull = cdfVector(xgrid, xall, fullcdf);
%toc



%%
%p.numBoots = 10^3;

cdfgridBootMat = nan(p.numBoots, numel(xgrid));

tic
for k = 1:p.numBoots
    ind = randsample(1:N, N, true); 
    tmpDO = DO(ind);
    [sortedUniqueX, cdfval] = poolingDO(tmpDO);
    xSorted = xgrid;
    cdfgridb = cdfVector(xSorted, sortedUniqueX, cdfval);
    cdfgridBootMat(k, :) = cdfgridb;
end
toc    

%%
ciu = quantile(cdfgridBootMat, 0.975);
cil = quantile(cdfgridBootMat, 0.025);

%figure;, 
%plot(xgrid, ciu)
%hold on 
%plot(xgrid, cil)

figout{4} = figure;
errbaru = ciu-cdfgridFull';
errbarl = cdfgridFull'-cil;
s1 = shadedErrorBarV2(xgrid, cdfgridFull, [errbaru; errbarl], 'lineprops', '-r');
s1.mainLine.Visible='off';

hold on

xx =  DO(1).uniqueX;
cdfval = DO(1).cdfval;
stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])
grid on
ax = gca; ax.YLim=[0, 1];

hold on
for k = 2:numel(DO)
    %[f, x] = ecdf(DASdistns{k}); plot(x, f)
    %ecdf(DASdistns{k})
    xx =  DO(k).uniqueX;
    cdfval = DO(k).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])

end

xlabel(varName)
ylabel('Cumulative frequency')
title([varName, '-', ctlName])

%legend('72M 1', '72M 2', '72M 3', 'Location', 'Southeast')
%legend('32M 1', '32M 2', '32M 3', 'Location', 'Southeast')

%
[xx, cdfval] = poolingDO(DO);
stairs([xx(1); xx], [0; cdfval], 'Color', 'r', 'LineWidth', 2)



%%  trt

DO = DO_trt;

% DO_ctrl
N = numel(DO);
[xall, fullcdf] = poolingDO(DO);

minx = min(xall); maxx = max(xall);
xgrid = linspace(minx, maxx, 1001);
numel(xgrid);

%tic
cdfgridFull = cdfVector(xgrid, xall, fullcdf);
%toc 


%%

cdfgridBootMat2 = nan(p.numBoots, numel(xgrid));

tic
for k = 1:p.numBoots
    ind = randsample(1:N, N, true); 
    tmpDO = DO(ind);
    [sortedUniqueX, cdfval] = poolingDO(tmpDO);
    xSorted = xgrid;
    cdfgridb = cdfVector(xSorted, sortedUniqueX, cdfval);
    cdfgridBootMat2(k, :) = cdfgridb;
end
toc    

%%
ciu2 = quantile(cdfgridBootMat2, 0.975);
cil2 = quantile(cdfgridBootMat2, 0.025);

%figure;, 
%plot(xgrid, ciu)
%hold on 
%plot(xgrid, cil)

figout{5} = figure;
errbaru = ciu2-cdfgridFull';
errbarl = cdfgridFull'-cil2;
s1 = shadedErrorBarV2(xgrid, cdfgridFull, [errbaru; errbarl], 'lineprops', '-r');
s1.mainLine.Visible='off';

hold on

xx =  DO(1).uniqueX;
cdfval = DO(1).cdfval;
stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])
grid on
ax = gca; ax.YLim=[0, 1];

hold on
for k = 2:numel(DO)
    %[f, x] = ecdf(DASdistns{k}); plot(x, f)
    %ecdf(DASdistns{k})
    xx =  DO(k).uniqueX;
    cdfval = DO(k).cdfval;
    stairs([xx(1); xx], [0; cdfval], 'Color', [.5 .5 .5])

end

xlabel(varName)
ylabel('Cumulative frequency')
title([varName, '-', trtName])

%legend('72M 1', '72M 2', '72M 3', 'Location', 'Southeast')
%legend('32M 1', '32M 2', '32M 3', 'Location', 'Southeast')

%
[xx, cdfval] = poolingDO(DO);
stairs([xx(1); xx], [0; cdfval], 'Color', 'r', 'LineWidth', 2)


%% together shadedbar

% DO_ctrl 
[xall, fullcdf] = poolingDO(DO_ctrl);
minx = min(xall); maxx = max(xall);
xgrid1 = linspace(minx, maxx, 1001);
numel(xgrid1);

cdfgridFull1 = cdfVector(xgrid1, xall, fullcdf);


% DO_trt 
[xall, fullcdf] = poolingDO(DO_trt);
minx = min(xall); maxx = max(xall);
xgrid2 = linspace(minx, maxx, 1001);
numel(xgrid2);

cdfgridFull2 = cdfVector(xgrid2, xall, fullcdf);


errbaru2 = ciu2-cdfgridFull2';
errbarl2 = cdfgridFull2'-cil2;


figout{6} = figure;
s1 = shadedErrorBarV2(xgrid1, cdfgridFull1, [errbaru; errbarl], 'lineprops', '-b');
s1.mainLine.LineWidth=2;

hold on
s2 = shadedErrorBarV2(xgrid2, cdfgridFull2, [errbaru2; errbarl2], 'lineprops', '-r');
s2.mainLine.LineWidth=2;
hold off
grid on
ax = gca;
ax.YLim = [0,1];

legend([s1.mainLine s2.mainLine], ctlName, trtName, 'Location', 'northeast')

xlabel(varName)
ylabel('Cumulative frequency')
title([ctlName, ' vs. ', trtName])


%%  saveas

if p.save

outFolder = p.outputDir;
if ~isdir(outFolder); mkdir(outFolder); end


saveas(figout{1}, fullfile(outFolder, 'fig1_ctrlCDFs.png'), 'png')
saveas(figout{2}, fullfile(outFolder, 'fig2_trtCDFs.png'), 'png')
saveas(figout{3}, fullfile(outFolder, 'fig3_ctrltrtCDFs.png'), 'png')
saveas(figout{4}, fullfile(outFolder, 'fig4_ctrlCI.png'), 'png')
saveas(figout{5}, fullfile(outFolder, 'fig5_trtCI.png'), 'png')
saveas(figout{6}, fullfile(outFolder, 'fig6_CIoverlayed.png'), 'png')


saveas(figout{1}, fullfile(outFolder, 'fig1_ctrlCDFs.fig'), 'fig')
saveas(figout{2}, fullfile(outFolder, 'fig2_trtCDFs.fig'), 'fig')
saveas(figout{3}, fullfile(outFolder, 'fig3_ctrltrtCDFs.fig'), 'fig')
saveas(figout{4}, fullfile(outFolder, 'fig4_ctrlCI.fig'), 'fig')
saveas(figout{5}, fullfile(outFolder, 'fig5_trtCI.fig'), 'fig')
saveas(figout{6}, fullfile(outFolder, 'fig6_CIoverlayed.fig'), 'fig')

end

%%
disp('=======')
disp('distnObj_plotCDFs_BootCI is done!')
if p.save
    disp('Outputs saved in:')
    disp(p.outputDir)
end

disp('====  close all')


end



%