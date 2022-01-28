function  newTracks=addSphericalCoord(obj,tracks)
   % Philippe Roudot 2017
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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
   for fIdx=1:length(obj)
       track=obj(fIdx);
       try
           track.addprop('azimuth');      
           track.addprop('elevation');    
           track.addprop('rho');          
       catch
       end
       [track.azimuth,track.elevation,track.rho]=cart2sph(det.xCoord(:,1), ... 
                                           det.yCoord(:,1), ...
                                           det.zCoord(:,1));
   
   end   
end
