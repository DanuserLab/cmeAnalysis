%[idx1, idx2] = colocalizationLAP(X1, X2, R) determines colocalization between two
% sets of points based on mutual proximity via bipartite graph matching.
%
% Inputs:
%      X1 : first set of points, up to 3D: [x1 y1 z1], NxDim format
%      X2 : second set of points
%
% Outputs:
%    idx1 : index of matched points in the first set 
%    idx2 : index of corresponding points in the second set
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

% Francois Aguet, 2013

function [idx1, idx2] = colocalizationLAP(X1, X2, R)

D = createSparseDistanceMatrix(X1, X2, R);
[link12, ~] = lap(D, [], [], 1);

n1 = size(X1,1);
n2 = size(X2,1);
link12 = link12(1:n1);
matchIdx = link12<=n2;
idx1 = find(matchIdx);
idx2 = double(link12(matchIdx));
