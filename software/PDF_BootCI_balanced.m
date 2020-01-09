function [fig, pdfgridBootMat, pts] = ...
    PDF_BootCI_balanced(cellData, numBoots, bw, title0, xlab, PosSuppArg)
% PDF_BootCI_balanced COMPUTE bootstrap confidence interval marginally
% using balanced pooling scheme.
%
% 2018/01/20, Jungsik Noh
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

%%
%xf = cell2mat(cellData);
xf = balancedPooling(cellData);
[ff,pts,bw0] = ksdensity(xf, 'NumPoints', 401, 'Bandwidth', bw, PosSuppArg{:});  % bw0 is the optimalBW for the original data.

N = numel(cellData);


%% bootstrap
%numBoots = 10^3;

pdfgridBootMat = nan(numBoots, numel(pts));

tic
parfor k = 1:numBoots
    ind = randsample(1:N, N, true); 
    tmpDO = cellData(ind);
    %bootPooledSample = cell2mat(tmpDO);
    bootPooledSample = balancedPooling(tmpDO);
    
    [f,~] = ksdensity(bootPooledSample, pts, 'BandWidth', bw0, PosSuppArg{:});
    pdfgridBootMat(k, :) = f;
end
toc    


%% plot
ciu = quantile(pdfgridBootMat, 0.975);
cil = quantile(pdfgridBootMat, 0.025);

fig = figure;

hold on
for k = 1:10
    y = pdfgridBootMat(k, :);
    %[f,xi] = ksdensity(x, pts, 'BandWidth', bw0);
    p1 = plot(pts, y, 'Color', [.5 .5 .5]);
end

errbaru = ciu-ff;
errbarl = ff-cil;
s1 = shadedErrorBarV2(pts, ff, [errbaru; errbarl], 'lineprops', '-r');
%s1.mainLine.Visible='off';
s1.mainLine.LineWidth = 2; 
%s1.patch.FaceAlpha = 0.3;


grid on

xlabel(xlab)
ylabel('Probability Density')

title([title0, ' bandWidth: ', num2str(round(bw0,2))])
legend(p1, 'Bootstrapped (10 among all)')
 

end
