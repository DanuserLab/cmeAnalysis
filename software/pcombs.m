function C = pcombs(v,keepDiag)
% This function generate all pairwise combinations of a vector v of length
% n. The number of possible pairwise combinations is n x (n - 1) / 2 with
% keepDiag == false (default), n x (n + 1) / 2 otherwise. The elements of v
% must be all different.
%
% e.g.:
% pcombs([1 2 4 3], false)
% ans = 
%   1  2
%   1  4
%   1  3
%   2  4
%   2  3
%   4  3
%
% pcombs([1 2 4 3], true)
% ans =
%   1  1
%   1  2
%   1  4
%   1  3
%   2  2
%   2  4
%   2  3
%   4  4
%   4  3
%   3  3
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

if nargin < 2 || isempty(keepDiag)
    keepDiag = false;
end

assert(islogical(keepDiag));

assert(length(v) == length(unique(v)));

[I J] = meshgrid(v, v');
k = ~keepDiag;
I = triu(I, k);
J = triu(J, k);
C = horzcat(nonzeros(J), nonzeros(I));

end