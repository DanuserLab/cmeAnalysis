function sem=SEM_calc(vect, CI)
% SEM_calc - standard error of the mean, confidence interval
%
% function sem=SEM_calc(vect, CI) 
%
% Purpose
% Calculate the standard error the mean to a given confidence
% interval (CI). Note that nans do not contribute to the
% calculation of the sample size and are ignored for the SD
% calculation. Output of this function has been checked against
% known working code written in R. 
%
% Inputs
% - vect: A vector upon which the SEM will be calculated. Note that
%         if vect is a matrix then we calculate one SEM for each
%         column. 
%
% - CI [optional]: a p value for a different 2-tailed interval. e.g. 0.01
%   This is a 2-tailed interval.
%
% Outputs
% sem - the standard error of the mean. So to plot the interval it's mu-sem
% to mu+sem. 
%
% Example - plot a 1% interval [rather than the default %5]
% r=randn(1,30);
% S=SEM_calc(r,0.01);
% hist(r)
% hold on
% plot(mean(r), mean(ylim),'r*')
% plot([mean(r)-S,mean(r)+S], [mean(ylim),mean(ylim)],'r-')
% hold off
%
% Rob Campbell 
%
% Also see - tInterval_Calc, norminv
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

error(nargchk(1,2,nargin))

if isvector(vect)
  vect=vect(:);
end


if nargin==1
  stdCI = 1.96 ; 
elseif nargin==2
  CI = CI/2 ; %Convert to 2-tail
  stdCI = abs(norminv(CI,0,1)) ;
end

sem = ( (nanstd(vect)) ./ sqrt(sum(~isnan(vect))) ) * stdCI ;    




