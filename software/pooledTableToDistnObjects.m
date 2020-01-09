function distnObjects = pooledTableToDistnObjects(pooledTable)
% pooledTableToDistnObjects converts pooledDataMatrix into a structure
% of the CDF objects.
%
% Input: 
%   - pooledDataMatrix has 2 columns: the 1st is for data points per
% experiment; the 2nd is for experiment ID.
%     Eg.  DAS      Movie_num
%           0.014       1
%           -10.5       1
%           ...         ...
%           2.5         12
%
%
% J Noh, 2018/01/20
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


%% nan check / subjectNames

pooledDataMatrix = table2array(pooledTable(:, 1:2));
if (size(pooledTable, 2) > 2)
    subjectNames = unique(pooledTable{:, 3});
else
    subjectNames = unique(pooledTable{:, 2});
end


disp('== num of nans in each column ==')
sum(isnan(pooledDataMatrix))


%% subjectID

%inunitDatapoints = pooledDataMatrix(:, 1);
subjectID = pooledDataMatrix(:, 2); 

disp('== subjectID ==')
tabulate(subjectID)
%tbl = tabulate(subjectID);

[subjectID_us, sorti] = sort(unique(subjectID));
disp('== unique & sorted subjectID ==')
disp(subjectID_us)

numExperiments = numel(subjectID_us);
disp(['== Number of subjects: ', num2str(numExperiments), ' =='])


%% indexing of experiment IDs

inds = cell(numExperiments, 1);

for i = 1:numExperiments
    inds{i} = (subjectID == subjectID_us(i));
end

%% make distnObjects

distnObjects = struct('subjectID', cell(numExperiments, 1));

for i=1:numExperiments
    inunitDatapoints = pooledDataMatrix(inds{i}, 1);
    xx = sort(inunitDatapoints);
    tmp = rle(xx);
    uniqueX = tmp(1:2:end);
    rlen = tmp(2:2:end);
    uniqueInd = cumsum(rlen);
    cumfreq = (1:numel(xx))./numel(xx);
    cdfval = cumfreq(uniqueInd);
    distnObjects(i).uniqueX = uniqueX(:);
    distnObjects(i).cdfval = cdfval(:);
    
    distnObjects(i).subjectID = subjectID_us(i);    
    distnObjects(i).subjectName = subjectNames(sorti(i));
end


%% output massage

nrc = size(pooledTable);
disp(['== A pooled Table of size ', num2str(nrc),  ' was converted to a struct with ', ...
    num2str(numExperiments), ' Distributional Objects.'])


end
