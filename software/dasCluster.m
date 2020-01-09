

function [idx, Cx] = dasCluster(DAS_all, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('correct13', false, @islogical);
ip.addParameter('correct13_GM', false, @islogical);
ip.parse(DAS_all, pm, varargin{:});


DAS_all = ip.Results.DAS_all;
correct13 = ip.Results.correct13;
num_condition = max(size(DAS_all));
correct13_GM = ip.Results.correct13_GM;
pm = ip.Results.pm;
%==========================================================================
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
idx = cell(num_condition,1); 
Cx = cell(num_condition,1); 
num_clus = 3;
%==========================================================================
disp('clustering ...')
for i_condition = 1:num_condition
    X = DAS_all{i_condition}.DAS_var;
    Y = DAS_all{i_condition}.DAS;
    Z = DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3;
    if i_condition == 1
        X1 = X; Y1 = Y; Z1 = Z;
        mu = [mean(X1,'omitnan'), mean(Y1,'omitnan'), mean(Z1,'omitnan')];
        sigma = [std(X1,'omitnan'), std(Y1,'omitnan'), std(Z1,'omitnan')];
    end
    
    if (i_condition == 1)
        Z_norm=[normalize(X,'zscore'),normalize(Y,'zscore'),normalize(Z,'zscore')];
        [idx{i_condition},Cx{i_condition}] = kmeans(Z_norm,num_clus,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
    elseif (i_condition > 1) && (pm.kmeans_ctrl_only == false)
        Z_norm=[normalize(X,'zscore'),normalize(Y,'zscore'),normalize(Z,'zscore')];
        [idx{i_condition},Cx{i_condition}] = kmeans(Z_norm,num_clus,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
    elseif (i_condition > 1) && (pm.kmeans_ctrl_only == true)
        dist_C = zeros(max(size(DAS_all{i_condition}.DAS)),num_clus);
        Z_norm=[(X-mu(1))/sigma(1),(Y-mu(2))/sigma(2),(Z-mu(3))/sigma(3)];
        for i_c = 1:num_clus
         dist_C(:,i_c) = sum((Z_norm-Cx{1}(i_c,:)).*(Z_norm-Cx{1}(i_c,:)),2);
        end
        [~,idx{i_condition}] = min(dist_C,[],2);
          Cx{i_condition} = Cx{1};
%           dist_C = zeros(max(size(DAS_all{i_condition}.DAS)),num_clus);
%         c_tem = Cx{1} .* sigma+mu;
%         Z_norm=[X,Y,Z];
%         for i_c = 1:num_clus
%          dist_C(:,i_c) = sum((Z_norm-c_tem(i_c,:)).*(Z_norm-c_tem(i_c,:)),2);
%         end
%         [~,idx{i_condition}] = min(dist_C,[],2);
%           Cx{i_condition} = Cx{1};
    end
    %[idx{i_condition},Cx{i_condition}] = kmeans(Z_norm,num_clus,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
    %----------------------------------------------------
     i_ccp = ones(num_clus,1)*2;
     LT_mean = zeros(num_clus,1);
     Imax_mean = zeros(num_clus,1);
     for i_c = 1:num_clus
      LT_mean(i_c) = mode(DAS_all{i_condition}.LT(idx{i_condition}==i_c));
      Imax_mean(i_c) = mean(DAS_all{i_condition}.MaxI(idx{i_condition}==i_c));
     end
     [~,i_LT_mean_max] = max(LT_mean);
     [~,i_Imax_mean_min] = min(Imax_mean);
      i_ccp(i_LT_mean_max) = 1;
      i_ccp(i_Imax_mean_min) = 3;
  %-----------------------------------------
  Cx_temp = Cx{i_condition};
  id_temp = cell(num_clus,1);
  for i_c = 1:num_clus
      id_temp{i_c} = idx{i_condition}==i_c;
  end
  for i_c = 1:num_clus     
      idx{i_condition}(id_temp{i_c}) = i_ccp(i_c);
      Cx{i_condition}(i_ccp(i_c),:) = Cx_temp(i_c,:);
  end
  %-----------------------------------------
  if correct13 == true
   if (i_condition == 1) || (pm.kmeans_ctrl_only == false)
   id_temp = (idx{i_condition} == 1) | (idx{i_condition} == 3);
   Z_norm=[normalize(X(id_temp),'zscore'),normalize(Y(id_temp),'zscore')];
   if correct13_GM == false
   [idx_temp,C_temp] = kmedoids(Z_norm,2,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
   else
   GMM = fitgmdist(Z_norm,2,'Options',statset('Display','final','MaxIter',1000),'Replicates',10,'Start','randSample');
   idx_temp = cluster(GMM,Z_norm);
   end
   if mean(Z_norm(idx_temp==1,2)) > mean(Z_norm(idx_temp==2,2))
       idx_temp(idx_temp == 2) = 3;
       idx{i_condition}(id_temp) = idx_temp;
   else
       idx_temp = idx_temp-1;
       idx_temp(idx_temp == 0) = 3;
       idx{i_condition}(id_temp) = idx_temp;
       C_temp_store = C_temp;
       C_temp(1,:) = C_temp_store(2,:);
       C_temp(2,:) = C_temp_store(1,:);
   end
   Cx{i_condition}(1,:) = C_temp(1,:);
   Cx{i_condition}(3,:) = C_temp(2,:);
   end
  end
  %-----------------------------------------
    fprintf('condition# %4.0f finished \n',i_condition)
end


    
    
