function Output=rle(Input)
L=length(Input);
j=1;
k=1;
i=1;
while i<2*L
    comp=1;
    for j=j:L
        if j==L 
            break
        end;  
         if Input(j)==Input(j+1)
            comp=comp+1;
        else
            break
        end;
    end;
        Output(k+1)=comp;
        Output(k)=Input(j);
        if j==L && Input(j-1)==Input(j) 
            break
        end;  
        i=i+1;
        k=k+2;
        j=j+1;
        if j==L 
            if mod(L,2)==0 
            Output(k+1)=1;
            Output(k)=Input(j);
            else
            Output(k+1)=1;    
            Output(k)=Input(j);
                       end;
             break
        end;
    end; 
% Function of RLE (Run Length Encoding) was used on program of gray level image compression .
% Authors : Said Bourezg  - Derbel Abd Elhak
% Electronics Engineer  option:communication .
% Date : 05.26.2009
% Filename rle.m (Matlab)
% This function is part of my undergraduate project in M'sila university, Algeria.
% Adress:                          Said BOUREZG
%                               Elbassatine street
%                                 28038 Tarmount
%                               M'sila --- Algeria 
% Email:  said.bourezg@yahoo.fr
% Mobile: +213 796 018049 
% If you can improve this code furtherly, please let me know. Thanks
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
