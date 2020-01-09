function [uniqueX, cdfval] = poolingDO(DOstruct)
% [uniqueX, cdfval] = poolingDO(DOstruct)
%
% Jungsik Noh, 02/20/2017
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

%% add relfreq

tmpDO = DOstruct;
len = numel(DOstruct);

for i=1:len
    cumfreq = DOstruct(i).cdfval;
    relfreq = [cumfreq(1); diff(cumfreq)];
    tmpDO(i).relfreq = relfreq(:);
end

xxcell = {tmpDO.uniqueX}';
relfreqcell = {tmpDO.relfreq}';
pxx = cell2mat(xxcell);
prelfreq = cell2mat(relfreqcell);

[xx, ord] = sort(pxx);
freq0 = prelfreq(ord);
freq1 = freq0./sum(freq0);
%freq1 = freq0./len;
tmp = rle(xx);
uniqueX = tmp(1:2:end);
rlen = tmp(2:2:end);
uniqueInd = cumsum(rlen);
cumfreq = cumsum(freq1);
cdfval = cumfreq(uniqueInd);

uniqueX = uniqueX(:);
cdfval = cdfval(:);


end
