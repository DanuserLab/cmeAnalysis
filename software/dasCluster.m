function [idx, Cx] = dasCluster(DAS_all, pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.parse(DAS_all, pm, varargin{:});

DAS_all = ip.Results.DAS_all;
num_condition = max(size(DAS_all));
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
if pm.dist_perct > 0
    Z_norm_save = cell(num_condition,1);
end
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
        [idx{i_condition},Cx{i_condition}] = kmedoids(Z_norm,num_clus,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
    elseif (i_condition > 1) 
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
    Z_norm_save{i_condition} = Z_norm;
    %[idx{i_condition},Cx{i_condition}] = kmeans(Z_norm,num_clus,'Start','sample','Options',statset('MaxIter',10000),'Replicates',10);
    %----------------------------------------------------
    %fprintf('condition# %4.0f finished \n',i_condition)
end
for i_condition = 1:num_condition
     if (i_condition == 1) 
     i_ccp = ones(num_clus,1)*1;
     %LT_mean = zeros(num_clus,1);
     Imax_mean = zeros(num_clus,1);
     DAS_mean = zeros(num_clus,1);
     for i_c = 1:num_clus
      %LT_mean(i_c) = median(DAS_all{i_condition}.LT(idx{i_condition}==i_c));
      Imax_mean(i_c) = median(DAS_all{i_condition}.MaxI(idx{i_condition}==i_c));
      DAS_mean(i_c) = median(DAS_all{i_condition}.DAS(idx{i_condition}==i_c));
     end
     [~,i_Imax_mean_min] = min(Imax_mean);
     [~,DAS_mean_max] = max(DAS_mean);
      i_ccp(i_Imax_mean_min) = 3;
      i_ccp(DAS_mean_max) = 2;
     end
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
end

if pm.correct13 == true
    for i_condition = 1:num_condition
        idx_tem = idx{i_condition};
        idx{i_condition}(idx_tem == 1) = 3;
        idx{i_condition}(idx_tem == 3) = 1;
        Cx_tem = Cx{i_condition};
        Cx{i_condition}(1,:) = Cx_tem(3,:);
        Cx{i_condition}(3,:) = Cx_tem(1,:);
    end
elseif pm.correct12 == true
    for i_condition = 1:num_condition
        idx_tem = idx{i_condition};
        idx{i_condition}(idx_tem == 1) = 2;
        idx{i_condition}(idx_tem == 2) = 1;
        Cx_tem = Cx{i_condition};
        Cx{i_condition}(1,:) = Cx_tem(2,:);
        Cx{i_condition}(2,:) = Cx_tem(1,:);
    end
elseif pm.correct23 == true
     for i_condition = 1:num_condition
        idx_tem = idx{i_condition};
        idx{i_condition}(idx_tem == 2) = 3;
        idx{i_condition}(idx_tem == 3) = 2;
        Cx_tem = Cx{i_condition};
        Cx{i_condition}(2,:) = Cx_tem(3,:);
        Cx{i_condition}(3,:) = Cx_tem(2,:);
    end
end
if pm.dist_perct > 0
    for i_condition = 1:num_condition
       dist_C = zeros(max(size(DAS_all{i_condition}.DAS)),num_clus);
        for i_c = 1:num_clus
         dist_C(:,i_c) = sum((Z_norm_save{i_condition}-Cx{i_condition}(i_c,:)).*(Z_norm_save{i_condition}-Cx{i_condition}(i_c,:)),2);
        end
    dist_C_diff = abs(dist_C(:,1) - dist_C(:,3));
    dist_C_diff_thre = prctile(dist_C_diff((idx{i_condition}==1)|(idx{i_condition}==3)),pm.dist_perct);
    id_um = false(size(idx{i_condition}));
    id_um(((idx{i_condition}==1)|(idx{i_condition}==3))&(dist_C_diff<dist_C_diff_thre)) = true;
    idx{i_condition}(id_um) = 0;
    end
end


    
    
