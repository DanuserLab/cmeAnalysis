function [cdfVec, indVec] = cdfVector(xSorted, sortedUniqueX, cdfval)
%
% Jungsik Noh, 02/25/2017
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

M = numel(xSorted);
indVec = nan(M, 1);
cdfVec = ones(M, 1);

imax = numel(sortedUniqueX);
i = 1;

for j = 1:M
    while (i < imax) && (sortedUniqueX(i) <= xSorted(j))
        i=i+1;
    end
    
    if (i == imax); break; end
    if (i == 1)
        indVec(j) = 0;
        cdfVec(j) = 0;
    else
        indVec(j) = i-1;
        cdfVec(j) = cdfval(i-1);
    end
end

%
end
