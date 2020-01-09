function D = createSparseDistanceMatrix(M,N,threshold,varargin)
% calculates and stores all distances below threshold into a sparse matrix
%
% SYNOPSIS   D=createSparseDistanceMatrix(M,N,threshold)
%            D=createSparseDistanceMatrix(M,N,threshold,eps)
%
% INPUT      M and N are the matrices containing the set of point coordinates.
%            They have to be doubles.
%            M and N can represent point positions in 1, 2 and 3D, as follows.
%            
%            In 1D: M=[ y1        and   N=[ y1
%                       y2                  y2
%                       ...                 ... 
%                       ym ]                yn ]
%
%                   Distances: dij = yj-yi
%
%            In 2D:
%                   M=[ y1 x1     and   N=[ y1 x1
%                       y2 x2               y2 x2
%                        ...                 ...
%                       ym xm ]             yn xn ]
%
%                   Distances: dij = sqrt( (yj-yi)^2 + (xj-xi)^2 )
%
%            In 3D:
%                   M=[ y1 x1 z1  and   N=[ y1 x1 z1
%                       y2 x2 z2            y2 x2 z2
%                         ...                ...
%                       ym xm zm ]          yn xn zn ]
%
%                   Distances: dij = sqrt( (yj-yi)^2 + (xj-xi)^2 + (zj-zi)^2 )
%
%
%            threshold : only the distances dij between the two set of points 
%                        M and N which are below the threshold are stored in 
%                        the (sparse) distance matrix
%
%            eps       : (optional, default value eps = 1e-10)
%                        By definition, sparse matrices contain only non-zero 
%                        elements.
%                        The sparse matrix returned by the function contains 
%                        only the distances dij <= threshold. For this reason,
%                        any distance dij > threshold, which is therefore NOT 
%                        stored in D, is considered to be equal to 0 by Matlab. 
%                        To be able to distinguish the two cases: dij = 0 from
%                        dij > threshold, all zero distances are replaced in D 
%                        by a small number, eps (default eps=1e-10).
%                        Thus the command:
% 
%                           find(D==eps) returns the indices of all zero distances,
%                                        whereas the command:
%                           find(D==0)   returns the indices of all distances > threshold.
% 
% OUTPUT      D        : sparse distance matrix. 
%
% Sebastien Besson, July 2011
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

% Input check
ip =inputParser;
ip.addRequired('M',@(x)isa(x,'double'));
ip.addRequired('N',@(x)isa(x,'double'));
ip.addRequired('threshold',@isscalar);
ip.addOptional('epsilon',1e-10,@isscalar);
ip.parse(M,N,threshold,varargin{:})

% handle NaNs - kdTree will "forget" to report a few distances if there are
% NaNs present in the input
finiteM = find(all(isfinite(M),2));
finiteN = find(all(isfinite(N),2));

% Query the points below the threshold using the KDTree
[points,distances] = KDTreeBallQuery(N(finiteN,:),M(finiteM,:),threshold);

% Generate the list of indices to create the sparse matrix
% points1=points;
% nzInd=find(~cellfun(@isempty,points));
% for i=nzInd', points1{i}(:)=i;end

% a much faster version of the above -- jonas, 10/2012
% Create a vector with ones and zeros so that we can use cumsum to create
% number of times we need to repeat a given entry, and use it to index
% "toRepeat"
nRepeats = cellfun(@numel, points);
toRepeat = find(nRepeats);
index = zeros(sum(nRepeats),1);
index([1;cumsum(nRepeats(toRepeat(1:end-1)))+1])=1;

% Create the sparse matrix. Index into finiteM, finiteN to account for
% NaN-rows that have been removed.
if ~isempty(toRepeat)
    D = sparse(finiteM(toRepeat(cumsum(index))),...
        finiteN(vertcat(points{:})),...
        max(vertcat(distances{:}),ip.Results.epsilon),size(M,1),size(N,1));
else
    D = sparse(size(M,1),size(N,1));
end
