function [figout] = distnObj_plotPDFs_BootCI(dat_ctrl, dat_trt, ctlName, trtName, varName, varargin)
% distnObj_plotPDFs_BootCI VISUALIZE multiple distributional data objects 
% through probability density functions (PDFs) with bootstrapped confidence
% intervals.
% 
% Usage:
%   distnObj_plotPDFs_BootCI(dat_ctrl, dat_trt, 'condition1', 'condition2', ...
%       'TelomereLength', 'PositiveSupport', true)
%
% Input: 
%   - dat_ctrl is a table object of pooled observations with 2 columns: 
%       the 1st is for data points per experiment; the 2nd is for experiment ID.
%     Eg.  DAS      Movie_num
%           0.014       1
%           -10.5       1
%           ...         ...
%           2.5         12
%
% Options:
%   - PositiveSupport: whether the observations should be always positive.
%   Default is false.
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
ip.addParameter('outputDir', fullfile(pwd, 'distnObj_PDFs_BootCI'));
ip.addParameter('numBoots', 400)
ip.addParameter('BoundedSupport', false);
ip.addParameter('CommonBandWidth', NaN);
ip.parse(varargin{:});
p = ip.Results;

if p.BoundedSupport == false
    PosSuppArg = {'Support', 'unbounded'};
else
    PosSuppArg = {'Support', p.BoundedSupport, 'BoundaryCorrection', 'reflection'};
end

 
%%

dat_ctrl.Properties.VariableNames = {varName, 'subjectID'};
dat_trt.Properties.VariableNames = {varName, 'subjectID'};


%% boxplot for ctrl

figout = cell(10,1);

figout{1} = figure;
boxplot(dat_ctrl{:, 1}, dat_ctrl.subjectID)
figure(figout{1});
xlabel('subjectID')
ylabel(varName)
title(ctlName)
h = refline([0, median(dat_ctrl{:, 1})]);
legend(h, 'Median of pooled data')

% boxplot for trt
figout{2} = figure;
boxplot(dat_trt{:, 1}, dat_trt.subjectID)
figure(figout{2});
xlabel('subjectID')
ylabel(varName)
title(trtName)
h = refline([0, median(dat_trt{:, 1})]);
legend(h, 'Median of pooled data')


%% cell format

N1 = numel(unique(dat_ctrl.subjectID));
N2 = numel(unique(dat_trt.subjectID));


datCell_ctrl = cell(N1, 1);
datCell_trt = cell(N2, 1);

idvec = sort(unique(dat_ctrl.subjectID));
for k = 1:N1
    g = idvec(k);
    ind0 = (dat_ctrl.subjectID == g);
    datCell_ctrl{k} = dat_ctrl{ind0, 1};
end

idvec = sort(unique(dat_trt.subjectID));
for k = 1:N2
    g = idvec(k);
    ind0 = (dat_trt.subjectID == g);
    datCell_trt{k} = dat_trt{ind0, 1};
end



%% pdf ctrl

%xf = dat_ctrl.var;
xf = balancedPooling(datCell_ctrl);

[ff,pts,bw0] = ksdensity(xf, 'NumPoints', 401, PosSuppArg{:});
%figure, plot(pts, ff)
idvec = sort(unique(dat_ctrl.subjectID));

figout{3} = figure;
x = datCell_ctrl{1};
[f,xi,bw] = ksdensity(x, pts, PosSuppArg{:});
plot(xi, f, 'LineWidth', 0.5)

hold on
for k = 2:N1
    x = datCell_ctrl{k};
    [f,xi,bw] = ksdensity(x, pts, PosSuppArg{:});
    plot(xi, f, 'LineWidth', 0.5)
end
%legend(mat2cell(idvec, N1, 1))

p1 = plot(pts, ff, 'r', 'LineWidth', 2);
%legend(p1, 'pooled')
xlabel(varName)
ylabel('Probability Density')
title([varName, '-', ctlName])

legend(num2str(idvec))

%% pdf trt

%xf = dat_trt.var;
xf = balancedPooling(datCell_trt);

[ff,pts,bw0] = ksdensity(xf, 'NumPoints', 401, PosSuppArg{:});
%figure, plot(pts, ff)
idvec = sort(unique(dat_trt.subjectID)); 

figout{4} = figure;
x = datCell_trt{1};
[f,xi,bw] = ksdensity(x, pts, PosSuppArg{:});
plot(xi, f)

hold on
for k = 2:N2
    x = datCell_trt{k};
    [f,xi,bw] = ksdensity(x, pts, PosSuppArg{:});
    plot(xi, f)
end
%legend(mat2cell(idvec, N1, 1))

p1 = plot(pts, ff, 'r', 'LineWidth', 2);
%legend(p1, 'pooled')
xlabel(varName)

ylabel('Probability Density')
title([varName, '-', trtName])

legend(num2str(idvec))


%%
%%  2. Bootstrap Confidence Interval of PDFs from movie-wise variation
%%

%cellData = datCell_ctrl; 
%p.numBoots = 10

%xf = cell2mat(datCell_ctrl);
xf_ctrl = balancedPooling(datCell_ctrl);
[~,~,bwctrl0] = ksdensity(xf_ctrl, 'NumPoints', 401, PosSuppArg{:});

if isnan(p.CommonBandWidth)
    disp('== Default bandWidth is optimal only when underlying distribution is close to a Gaussian distribution.')
    bwctrl = bwctrl0;
else
    bwctrl = p.CommonBandWidth;
    bwctrl0 = p.CommonBandWidth;
end

disp(['bwctrl = ', num2str(bwctrl)])

[fctrl, pdfgridBootMat_ctrl, pts_ctrl] = ...
    PDF_BootCI_balanced(datCell_ctrl, p.numBoots, bwctrl, ctlName, varName, PosSuppArg);

[ftrt, pdfgridBootMat_trt, pts_trt] = ...
    PDF_BootCI_balanced(datCell_trt, p.numBoots, bwctrl, trtName, varName, PosSuppArg);

%
figout{5} = fctrl;
figout{6} = ftrt;


%% together
figout{7} = figure;

%xctrl = balancedPooling(datCell_ctrl);  
% xf_ctrl
[fc,ptsc,bw0] = ksdensity(xf_ctrl, pts_ctrl, 'Bandwidth', bwctrl, PosSuppArg{:});

ciuc = quantile(pdfgridBootMat_ctrl, 0.975);
cilc = quantile(pdfgridBootMat_ctrl, 0.025);
errbaru1 = ciuc-fc;
errbarl1 = fc-cilc;

%xtrt = cell2mat(datCell_trt);
xtrt = balancedPooling(datCell_trt);
[ft,ptst,bw0] = ksdensity(xtrt, pts_trt, 'Bandwidth', bwctrl, PosSuppArg{:});

ciut = quantile(pdfgridBootMat_trt, 0.975);
cilt = quantile(pdfgridBootMat_trt, 0.025);
errbaru2 = ciut-ft;
errbarl2 = ft-cilt;

grid on
hold on

s1 = shadedErrorBarV2(ptsc, fc, [errbaru1; errbarl1], 'lineprops', '-b');
s1.mainLine.LineWidth = 2; 

s2 = shadedErrorBarV2(ptst, ft, [errbaru2; errbarl2], 'lineprops', '-r');
s2.mainLine.LineWidth = 2; 


legend([s1.mainLine s2.mainLine], ctlName, trtName, 'Location', 'northeast')

xlabel(varName)
ylabel('Probability Density')
title([ctlName, ' vs. ', trtName, '. bandWidth: ', num2str(round(bwctrl,2))])


%%mobID
%% bw = bw0* 0.5

%cellData = datCell_ctrl; 
%p.numBoots = 10
% bwctrl = bwctrl * 0.5

%xf = cell2mat(datCell_ctrl);
%xf_ctrl = balancedPooling(datCell_ctrl);
%[~,~,bwctrl] = ksdensity(xf_ctrl, 'NumPoints', 401, PosSuppArg{:});

bwctrl = bwctrl0 * 0.5

[fctrl, pdfgridBootMat_ctrl, pts_ctrl] = ...
    PDF_BootCI_balanced(datCell_ctrl, p.numBoots, bwctrl, ctlName, varName, PosSuppArg);

[ftrt, pdfgridBootMat_trt, pts_trt] = ...
    PDF_BootCI_balanced(datCell_trt, p.numBoots, bwctrl, trtName, varName, PosSuppArg);

%
figout{8} = fctrl;
figout{9} = ftrt;


%% together
figout{10} = figure;

%xctrl = balancedPooling(datCell_ctrl);  
% xf_ctrl
[fc,ptsc,bw0] = ksdensity(xf_ctrl, pts_ctrl, 'Bandwidth', bwctrl, PosSuppArg{:});

ciuc = quantile(pdfgridBootMat_ctrl, 0.975);
cilc = quantile(pdfgridBootMat_ctrl, 0.025);
errbaru1 = ciuc-fc;
errbarl1 = fc-cilc;

%xtrt = cell2mat(datCell_trt);
xtrt = balancedPooling(datCell_trt);
[ft,ptst,bw0] = ksdensity(xtrt, pts_trt, 'Bandwidth', bwctrl, PosSuppArg{:});

ciut = quantile(pdfgridBootMat_trt, 0.975);
cilt = quantile(pdfgridBootMat_trt, 0.025);
errbaru2 = ciut-ft;
errbarl2 = ft-cilt;

grid on
hold on

s1 = shadedErrorBarV2(ptsc, fc, [errbaru1; errbarl1], 'lineprops', '-b');
s1.mainLine.LineWidth = 2; 

s2 = shadedErrorBarV2(ptst, ft, [errbaru2; errbarl2], 'lineprops', '-r');
s2.mainLine.LineWidth = 2; 


legend([s1.mainLine s2.mainLine], ctlName, trtName, 'Location', 'northeast')

xlabel(varName)
ylabel('Probability Density')
title([ctlName, ' vs. ', trtName, '. bandWidth: ', num2str(round(bwctrl,2))])



%% bw = bw0*2

%cellData = datCell_ctrl; 
%p.numBoots = 10   bwctrl = bwctrl * 2


%xf = cell2mat(datCell_ctrl);
%xf_ctrl = balancedPooling(datCell_ctrl);
%[~,~,bwctrl] = ksdensity(xf_ctrl, 'NumPoints', 401, PosSuppArg{:});

bwctrl = bwctrl0 * 2

[fctrl, pdfgridBootMat_ctrl, pts_ctrl] = ...
    PDF_BootCI_balanced(datCell_ctrl, p.numBoots, bwctrl, ctlName, varName, PosSuppArg);

[ftrt, pdfgridBootMat_trt, pts_trt] = ...
    PDF_BootCI_balanced(datCell_trt, p.numBoots, bwctrl, trtName, varName, PosSuppArg);

%
figout{11} = fctrl;
figout{12} = ftrt;


%% together
figout{13} = figure;

%xctrl = balancedPooling(datCell_ctrl);  
% xf_ctrl
[fc,ptsc,bw0] = ksdensity(xf_ctrl, pts_ctrl, 'Bandwidth', bwctrl, PosSuppArg{:});

ciuc = quantile(pdfgridBootMat_ctrl, 0.975);
cilc = quantile(pdfgridBootMat_ctrl, 0.025);
errbaru1 = ciuc-fc;
errbarl1 = fc-cilc;

%xtrt = cell2mat(datCell_trt);
xtrt = balancedPooling(datCell_trt);
[ft,ptst,bw0] = ksdensity(xtrt, pts_trt, 'Bandwidth', bwctrl, PosSuppArg{:});

ciut = quantile(pdfgridBootMat_trt, 0.975);
cilt = quantile(pdfgridBootMat_trt, 0.025);
errbaru2 = ciut-ft;
errbarl2 = ft-cilt;

grid on
hold on

s1 = shadedErrorBarV2(ptsc, fc, [errbaru1; errbarl1], 'lineprops', '-b');
s1.mainLine.LineWidth = 2; 

s2 = shadedErrorBarV2(ptst, ft, [errbaru2; errbarl2], 'lineprops', '-r');
s2.mainLine.LineWidth = 2; 


legend([s1.mainLine s2.mainLine], ctlName, trtName, 'Location', 'northeast')

xlabel(varName)
ylabel('Probability Density')
title([ctlName, ' vs. ', trtName, '. bandWidth: ', num2str(round(bwctrl,2))])


%% saveas

if p.save

outFolder = p.outputDir;
if ~isdir(outFolder); mkdir(outFolder); end


saveas(figout{1}, fullfile(outFolder, 'BP_ctrl.png'), 'png')
saveas(figout{2}, fullfile(outFolder, 'BP_trt.png'), 'png')
saveas(figout{3}, fullfile(outFolder, 'pdf_ctrl.png'), 'png')
saveas(figout{4}, fullfile(outFolder, 'pdf_trt.png'), 'png')
saveas(figout{5}, fullfile(outFolder, 'bootCI_ctrl_bw0.png'), 'png')
saveas(figout{6}, fullfile(outFolder, 'bootCI_trt_bw0.png'), 'png')
saveas(figout{7}, fullfile(outFolder, 'CI_together_bw0.png'), 'png')
saveas(figout{8}, fullfile(outFolder, 'bootCI_ctrl_bw0_half.png'), 'png')
saveas(figout{9}, fullfile(outFolder, 'bootCI_trt_bw0_half.png'), 'png')
saveas(figout{10}, fullfile(outFolder, 'CI_together_bw0_half.png'), 'png')
saveas(figout{11}, fullfile(outFolder, 'bootCI_ctrl_bw0_double.png'), 'png')
saveas(figout{12}, fullfile(outFolder, 'bootCI_trt_bw0_double.png'), 'png')
saveas(figout{13}, fullfile(outFolder, 'CI_together_bw0_double.png'), 'png')



saveas(figout{1}, fullfile(outFolder, 'BP_ctrl.fig'), 'fig')
saveas(figout{2}, fullfile(outFolder, 'BP_trt.fig'), 'fig')
saveas(figout{3}, fullfile(outFolder, 'pdf_ctrl.fig'), 'fig')
saveas(figout{4}, fullfile(outFolder, 'pdf_trt.fig'), 'fig')
saveas(figout{5}, fullfile(outFolder, 'bootCI_ctrl_bw0.fig'), 'fig')
saveas(figout{6}, fullfile(outFolder, 'bootCI_trt_bw0.fig'), 'fig')
saveas(figout{7}, fullfile(outFolder, 'CI_together_bw0.fig'), 'fig')
saveas(figout{8}, fullfile(outFolder, 'bootCI_ctrl_bw0_half.fig'), 'fig')
saveas(figout{9}, fullfile(outFolder, 'bootCI_trt_bw0_half.fig'), 'fig')
saveas(figout{10}, fullfile(outFolder, 'CI_together_bw0_half.fig'), 'fig')
saveas(figout{11}, fullfile(outFolder, 'bootCI_ctrl_bw0_double.fig'), 'fig')
saveas(figout{12}, fullfile(outFolder, 'bootCI_trt_bw0_double.fig'), 'fig')
saveas(figout{13}, fullfile(outFolder, 'CI_together_bw0_double.fig'), 'fig')

end

%%
disp('=======')
disp(['distnObj_plotPDFs_BootCI is done!'])
if p.save
    disp('Outputs saved in:')
    disp(p.outputDir)
end

disp('====  close all')

end


