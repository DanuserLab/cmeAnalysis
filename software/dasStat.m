
function [DAS_stat, idx] = dasStat(dir_alt, pm, DAS_all, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('DAS_all', @(x) iscell(x));
ip.addParameter('d_bw', 0.05, @isnumeric);
ip.addParameter('overwriteStat', false, @islogical);
ip.parse(dir_alt, pm, DAS_all, varargin{:});

%--------------------------------------------------------------------------
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
DAS_all = ip.Results.DAS_all;
dir_alt = ip.Results.dir_alt;
overwriteStat = ip.Results.overwriteStat;
pm = ip.Results.pm;
%--------------------------------------------------------------------------
num_condition = max(size(DAS_all));
DAS_stat = cell(num_condition,1);
num_clus = 3;
%--------------------------------------------------------------------------------------------------------
if (overwriteStat == true) || (~(exist(fullfile(dir_alt, 'DAS_stat_idx.mat'), 'file') == 2))
%================================================================================
for i_condition = 1:num_condition
DAS_stat{i_condition} = struct(...
    'd',[], 'pdf_d',[], 'dv', [], 'pdf_dv',[],'d2', [], 'pdf_d2',[],'d3', [], 'pdf_d3',[],'d32', [], 'pdf_d32',[],...
    'LT',[], 'pdf_LT',[],'Imax',[], 'pdf_Imax',[],...
    'LT1',[], 'pdf_LT1',[],'Imax1',[], 'pdf_Imax1',[],...
    'LT2',[], 'pdf_LT2',[],'Imax2',[], 'pdf_Imax2',[],...
    'LT3',[], 'pdf_LT3',[],'Imax3',[], 'pdf_Imax3',[],...
    'z_d_dv',[],'edges_d_dv',[],...
    'z_d_d32',[],'edges_d_d32',[],...
    'z_d32_dv',[],'edges_d32_dv',[],...
    'd_mod',[],'dv_mod',[],'d32_mod',[],...
    'prop_c',[], ...
    'LT_med',[], ...
    'Imax_med',[], ...
    'p_LT_med',[], ...
    'p_Imax_med',[], ...
    'ps_LT_med',[], ...
    'ps_Imax_med',[], ...
    'p_LT_mean',[], ...
    'p_Imax_mean',[], ...
    'ps_LT_mean',[], ...
    'ps_Imax_mean',[], ...
    'tot_init',[] ...
     );
end
[idx, Cx] = dasCluster(DAS_all,pm);
%================================================================================
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.DAS,'Normalization','pdf');
DAS_stat{i_condition}.d = h.BinEdges;DAS_stat{i_condition}.d(end) = [];DAS_stat{i_condition}.pdf_d = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.DAS_var,'Normalization','pdf');
DAS_stat{i_condition}.dv = h.BinEdges;DAS_stat{i_condition}.dv(end) = [];DAS_stat{i_condition}.pdf_dv = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.DAS_2,'Normalization','pdf');
DAS_stat{i_condition}.d2 = h.BinEdges;DAS_stat{i_condition}.d2(end) = [];DAS_stat{i_condition}.pdf_d2 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.DAS_3,'Normalization','pdf');
DAS_stat{i_condition}.d3 = h.BinEdges;DAS_stat{i_condition}.d3(end) = [];DAS_stat{i_condition}.pdf_d3 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3,'Normalization','pdf');
DAS_stat{i_condition}.d32 = h.BinEdges;DAS_stat{i_condition}.d32(end) = [];DAS_stat{i_condition}.pdf_d32 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
BinWidth_LT_all = cell(num_condition,1);
BinLim_LT_all = cell(num_condition,1);
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.LT,'Normalization','pdf');
BinWidth_LT_all{i_condition} = h.BinWidth;
BinLim_LT_all{i_condition} = [5, max(DAS_all{i_condition}.LT)];
h = histogram(DAS_all{i_condition}.LT,'Normalization','pdf','BinLimits',BinLim_LT_all{i_condition},'BinWidth',BinWidth_LT_all{i_condition});
DAS_stat{i_condition}.LT = h.BinEdges;DAS_stat{i_condition}.LT(end) = [];DAS_stat{i_condition}.pdf_LT = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
BinWidth_Imax_all = cell(num_condition,1);
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.MaxI,'Normalization','pdf');
DAS_stat{i_condition}.Imax = h.BinEdges;DAS_stat{i_condition}.Imax(end) = [];DAS_stat{i_condition}.pdf_Imax = h.Values;
BinWidth_Imax_all{i_condition} = h.BinWidth;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.LT(idx{i_condition}==1),'Normalization','pdf', 'BinWidth',BinWidth_LT_all{i_condition},'BinLimits',BinLim_LT_all{i_condition});
DAS_stat{i_condition}.LT1 = h.BinEdges;DAS_stat{i_condition}.LT1(end) = [];DAS_stat{i_condition}.pdf_LT1 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.MaxI(idx{i_condition}==1),'Normalization','pdf', 'BinWidth',BinWidth_Imax_all{i_condition} );
DAS_stat{i_condition}.Imax1 = h.BinEdges;DAS_stat{i_condition}.Imax1(end) = [];DAS_stat{i_condition}.pdf_Imax1 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.LT(idx{i_condition}==2),'Normalization','pdf' , 'BinWidth',BinWidth_LT_all{i_condition},'BinLimits',BinLim_LT_all{i_condition});
DAS_stat{i_condition}.LT2 = h.BinEdges;DAS_stat{i_condition}.LT2(end) = [];DAS_stat{i_condition}.pdf_LT2 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.MaxI(idx{i_condition}==2),'Normalization','pdf', 'BinWidth',BinWidth_Imax_all{i_condition} );
DAS_stat{i_condition}.Imax2 = h.BinEdges;DAS_stat{i_condition}.Imax2(end) = [];DAS_stat{i_condition}.pdf_Imax2 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.LT(idx{i_condition}==3),'Normalization','pdf', 'BinWidth',BinWidth_LT_all{i_condition},'BinLimits',BinLim_LT_all{i_condition});
DAS_stat{i_condition}.LT3 = h.BinEdges;DAS_stat{i_condition}.LT3(end) = [];DAS_stat{i_condition}.pdf_LT3 = h.Values;
hold off
close(gcf)
end
%--------------------------------------------------------------------------
for i_condition = 1:num_condition
figure;hold on;
h = histogram(DAS_all{i_condition}.MaxI(idx{i_condition}==3),'Normalization','pdf', 'BinWidth',BinWidth_Imax_all{i_condition} );
DAS_stat{i_condition}.Imax3 = h.BinEdges;DAS_stat{i_condition}.Imax3(end) = [];DAS_stat{i_condition}.pdf_Imax3 = h.Values;
hold off
close(gcf)
end
%================================================================================
range_1 = (-5.5:0.5:0);
range_2 = (-1.6:0.2:1); 
range_3 = (-3:0.5:3);
for i_condition = 1:num_condition
   Z=[DAS_all{i_condition}.DAS_var,DAS_all{i_condition}.DAS,DAS_all{i_condition}.DAS_3./sqrt(DAS_all{i_condition}.DAS_2).^3];
   %--------------------------------------------------------------------------
   DAS_stat{i_condition}.edges_d_dv = cell(2,1); 
       DAS_stat{i_condition}.edges_d_dv{2} = range_2;
       DAS_stat{i_condition}.edges_d_dv{1} = range_1;
   DAS_stat{i_condition}.z_d_dv = hist3(Z(:,1:2),'Edges',DAS_stat{i_condition}.edges_d_dv);
   %--------------------------------------------------------------------------
   DAS_stat{i_condition}.edges_d_d32 = cell(2,1); 
       DAS_stat{i_condition}.edges_d_d32{2} = range_2; 
       DAS_stat{i_condition}.edges_d_d32{1} = range_3;
   DAS_stat{i_condition}.z_d_d32 = hist3([Z(:,3),Z(:,2)],'Edges',DAS_stat{i_condition}.edges_d_d32);
   %--------------------------------------------------------------------------
   DAS_stat{i_condition}.edges_d32_dv = cell(2,1); 
       DAS_stat{i_condition}.edges_d32_dv{2} = range_3;
       DAS_stat{i_condition}.edges_d32_dv{1} = range_1;
   DAS_stat{i_condition}.z_d32_dv = hist3([Z(:,1),Z(:,3)],'Edges',DAS_stat{i_condition}.edges_d32_dv);
   %--------------------------------------------------------------------------
end
%================================================================================
for i_condition = 1:num_condition
   DAS_stat{i_condition}.d_mod = zeros(num_clus,1);
   DAS_stat{i_condition}.dv_mod = zeros(num_clus,1);
   DAS_stat{i_condition}.d32_mod = zeros(num_clus,1);
   for i_c = 1:num_clus
       id_temp = idx{i_condition}==i_c;
   Z=[DAS_all{i_condition}.DAS_var(id_temp),DAS_all{i_condition}.DAS(id_temp),DAS_all{i_condition}.DAS_3(id_temp)./sqrt(DAS_all{i_condition}.DAS_2(id_temp)).^3];
   %--------------------------------------------------------------------------
   z_temp = hist3(Z(:,1:2),'Edges',DAS_stat{i_condition}.edges_d_dv);
   [~,I] = max(z_temp(:));
   [I_row, I_col] = ind2sub(size(DAS_stat{i_condition}.z_d_dv),I);
   if I_row > max(size(DAS_stat{i_condition}.edges_d_dv{1}))
       I_row = max(size(DAS_stat{i_condition}.edges_d_dv{1}));
   end
   if I_col > max(size(DAS_stat{i_condition}.edges_d_dv{2}))
       I_col = max(size(DAS_stat{i_condition}.edges_d_dv{2}));
   end
   DAS_stat{i_condition}.d_mod(i_c) = DAS_stat{i_condition}.edges_d_dv{2}(I_col);
   DAS_stat{i_condition}.dv_mod(i_c) = DAS_stat{i_condition}.edges_d_dv{1}(I_row);
   %--------------------------------------------------------------------------
   z_temp = hist3([Z(:,3),Z(:,2)],'Edges',DAS_stat{i_condition}.edges_d_d32);
   [~,I] = max(z_temp(:));
   [I_row, I_col] = ind2sub(size(DAS_stat{i_condition}.z_d_dv),I);
   if I_row > max(size(DAS_stat{i_condition}.edges_d_d32{1}))
       I_row = max(size(DAS_stat{i_condition}.edges_d_d32{1}));
   end
   
   DAS_stat{i_condition}.d32_mod(i_c) = DAS_stat{i_condition}.edges_d_d32{1}(I_row);
   %--------------------------------------------------------------------------
   end
end
%================================================================================
for i_condition = 1:num_condition
DAS_stat{i_condition}.prop_c = zeros(3,1);
for i_c = 1:3
    DAS_stat{i_condition}.prop_c(i_c) = numel(idx{i_condition}(idx{i_condition}==i_c)) / numel(idx{i_condition});
end
num_movie = size(DAS_all{i_condition}.cellAreaTime,1);
DAS_stat{i_condition}.tot_init = zeros(num_movie,1);
for i_mov = 1:num_movie
    DAS_stat{i_condition}.tot_init(i_mov) = (numel(DAS_all{i_condition}.LT(DAS_all{i_condition}.MovieNum == i_mov))+DAS_all{i_condition}.Ib_num(i_mov)+DAS_all{i_condition}.Id_num(i_mov))/DAS_all{i_condition}.cellAreaTime(i_mov);
end                                           
end
%================================================================================
for i_condition = 1:num_condition
    DAS_stat{i_condition}.ps_LT_med = cell(3,1);
    DAS_stat{i_condition}.ps_Imax_med = cell(3,1);
    num_movie = size(DAS_all{i_condition}.cellAreaTime,1);
    x_median = zeros(num_movie,3); 
    y_median = zeros(num_movie,3);
    DAS_stat{i_condition}.LT_med = zeros(3,1);
    DAS_stat{i_condition}.Imax_med = zeros(3,1);
    if i_condition == 1
        DAS_stat{i_condition}.p_LT_med = [];
        DAS_stat{i_condition}.p_Imax_med = [];
    else
        DAS_stat{i_condition}.p_LT_med = zeros(3,1);
        DAS_stat{i_condition}.p_Imax_med = zeros(3,1);
    end
    for i_mov=1:num_movie        
           for i_c = 1:3 
            x=DAS_all{i_condition}.LT((idx{i_condition}==i_c)&(DAS_all{i_condition}.MovieNum==i_mov));
            x_median(i_mov,i_c)=median(x);
            y=DAS_all{i_condition}.MaxI((idx{i_condition}==i_c)&(DAS_all{i_condition}.MovieNum==i_mov));
            y_median(i_mov,i_c)=median(y);
           end
    end
    if i_condition == 1
        x_median_wt = x_median; y_median_wt = y_median;
    else
        for i_c = 1:3 
            [~,DAS_stat{i_condition}.p_LT_med(i_c)] = kstest2(x_median_wt(:,i_c),x_median(:,i_c));
            [~,DAS_stat{i_condition}.p_Imax_med(i_c)] = kstest2(y_median_wt(:,i_c),y_median(:,i_c));
        end
    end
    for i_c = 1:3 
        DAS_stat{i_condition}.LT_med(i_c) = mean(x_median(:,i_c),1);
        DAS_stat{i_condition}.Imax_med(i_c) = mean(y_median(:,i_c),1);
    end
%---------------------------------------
if i_condition > 1
for i_c = 1:3
    pos_neg = mean(x_median(:,i_c),1) > mean(x_median_wt(:,i_c),1);
    if DAS_stat{i_condition}.p_LT_med(i_c) > 0.05
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\rightarrow n.s.';
    elseif (DAS_stat{i_condition}.p_LT_med(i_c) <= 0.05) && (DAS_stat{i_condition}.p_LT_med(i_c) > 0.01)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\uparrow *';
       else
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\downarrow *';
       end
    elseif (DAS_stat{i_condition}.p_LT_med(i_c) <= 0.01) && (DAS_stat{i_condition}.p_LT_med(i_c) > 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\uparrow **';
       else
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\downarrow **';
       end
    elseif (DAS_stat{i_condition}.p_LT_med(i_c) <= 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\uparrow ***';
       else
           DAS_stat{i_condition}.ps_LT_med{i_c} = '\downarrow ***';
       end
    end

    pos_neg = mean(y_median(:,i_c),1) > mean(y_median_wt(:,i_c),1);
    if DAS_stat{i_condition}.p_Imax_med(i_c) > 0.05
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\rightarrow n.s.';
    elseif (DAS_stat{i_condition}.p_Imax_med(i_c) <= 0.05) && (DAS_stat{i_condition}.p_Imax_med(i_c) > 0.01)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\uparrow *';
       else
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\downarrow *';
       end
    elseif (DAS_stat{i_condition}.p_Imax_med(i_c) <= 0.01) && (DAS_stat{i_condition}.p_Imax_med(i_c) > 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\uparrow **';
       else
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\downarrow **';
       end
    elseif (DAS_stat{i_condition}.p_Imax_med(i_c) <= 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\uparrow ***';
       else
           DAS_stat{i_condition}.ps_Imax_med{i_c} = '\downarrow ***';
       end
    end
    
end
end
%---------------------------------------
end
%================================================================================
for i_condition = 1:num_condition
    DAS_stat{i_condition}.ps_LT_mean = cell(3,1);
    DAS_stat{i_condition}.ps_Imax_mean = cell(3,1);
    num_movie = size(DAS_all{i_condition}.cellAreaTime,1);
    x_mean = zeros(num_movie,3); 
    y_mean = zeros(num_movie,3);
    if i_condition == 1
        DAS_stat{i_condition}.p_LT_mean = [];
        DAS_stat{i_condition}.p_Imax_mean = [];
    else
        DAS_stat{i_condition}.p_LT_mean = zeros(3,1);
        DAS_stat{i_condition}.p_Imax_mean = zeros(3,1);
    end
    for i_mov=1:num_movie        
           for i_c = 1:3 
            x=DAS_all{i_condition}.LT((idx{i_condition}==i_c)&(DAS_all{i_condition}.MovieNum==i_mov));
            x_mean(i_mov,i_c)=mean(x);
            y=DAS_all{i_condition}.MaxI((idx{i_condition}==i_c)&(DAS_all{i_condition}.MovieNum==i_mov));
            y_mean(i_mov,i_c)=mean(y);
           end
    end
    if i_condition == 1
        x_mean_wt = x_mean; y_mean_wt = y_mean;
    else
        for i_c = 1:3 
            [~,DAS_stat{i_condition}.p_LT_mean(i_c)] = ranksum(x_mean_wt(:,i_c),x_mean(:,i_c));%kstest2
            [~,DAS_stat{i_condition}.p_Imax_mean(i_c)] = ranksum(y_mean_wt(:,i_c),y_mean(:,i_c));
        end
    end
%---------------------------------------
if i_condition > 1
for i_c = 1:3
    pos_neg = mean(x_mean(:,i_c),1) > mean(x_mean_wt(:,i_c),1);
    if DAS_stat{i_condition}.p_LT_mean(i_c) > 0.05
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\rightarrow n.s.';
    elseif (DAS_stat{i_condition}.p_LT_mean(i_c) <= 0.05) && (DAS_stat{i_condition}.p_LT_mean(i_c) > 0.01)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\uparrow *';
       else
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\downarrow *';
       end
    elseif (DAS_stat{i_condition}.p_LT_mean(i_c) <= 0.01) && (DAS_stat{i_condition}.p_LT_mean(i_c) > 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\uparrow **';
       else
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\downarrow **';
       end
    elseif (DAS_stat{i_condition}.p_LT_mean(i_c) <= 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\uparrow ***';
       else
           DAS_stat{i_condition}.ps_LT_mean{i_c} = '\downarrow ***';
       end
    end

    pos_neg = mean(y_mean(:,i_c),1) > mean(y_mean_wt(:,i_c),1);
    if DAS_stat{i_condition}.p_Imax_mean(i_c) > 0.05
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\rightarrow n.s.';
    elseif (DAS_stat{i_condition}.p_Imax_mean(i_c) <= 0.05) && (DAS_stat{i_condition}.p_Imax_mean(i_c) > 0.01)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\uparrow *';
       else
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\downarrow *';
       end
    elseif (DAS_stat{i_condition}.p_Imax_mean(i_c) <= 0.01) && (DAS_stat{i_condition}.p_Imax_mean(i_c) > 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\uparrow **';
       else
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\downarrow **';
       end
    elseif (DAS_stat{i_condition}.p_Imax_mean(i_c) <= 0.001)
       if pos_neg
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\uparrow ***';
       else
           DAS_stat{i_condition}.ps_Imax_mean{i_c} = '\downarrow ***';
       end
    end
    
end
end
%---------------------------------------
end
%================================================================================
n_bar=pm.n_bar;
test_target = cell(n_bar,1);
    for i_bar = 1:n_bar
        test_target{i_bar} = cell(num_condition,1);
    end
num_movie = zeros(num_condition,1);
for i_condition = 1:num_condition
        num_movie(i_condition) = size(DAS_all{i_condition}.cellAreaTime,1);
        for i_bar = 1:n_bar
        test_target{i_bar}{i_condition} = zeros(num_movie(i_condition),1);
        end
end
for i_condition = 1:num_condition
    for i_mov = 1:num_movie(i_condition)
            id_tem1 = (idx{i_condition}==1) & (DAS_all{i_condition}.MovieNum == i_mov);
            id_tem2 = (idx{i_condition}==2) & (DAS_all{i_condition}.MovieNum == i_mov);
            id_tem3 = (idx{i_condition}==3) & (DAS_all{i_condition}.MovieNum == i_mov);
            id_all = DAS_all{i_condition}.MovieNum == i_mov;
            
            test_target{1}{i_condition}(i_mov) = numel(DAS_all{i_condition}.LT(id_tem1))/DAS_all{i_condition}.cellAreaTime(i_mov);
            test_target{2}{i_condition}(i_mov) = numel(DAS_all{i_condition}.LT(id_tem2))/DAS_all{i_condition}.cellAreaTime(i_mov);
            test_target{3}{i_condition}(i_mov) = numel(DAS_all{i_condition}.LT(id_tem3))/DAS_all{i_condition}.cellAreaTime(i_mov);
            
            test_target{4}{i_condition}(i_mov) = (numel(DAS_all{i_condition}.LT(id_all))+DAS_all{i_condition}.Ib_num(i_mov))/DAS_all{i_condition}.cellAreaTime(i_mov);
            test_target{5}{i_condition}(i_mov) = numel(DAS_all{i_condition}.LT(id_tem1))/(numel(DAS_all{i_condition}.LT(id_tem1))+numel(DAS_all{i_condition}.LT(id_tem3)));
            
            if pm.plot_visitor == true
            test_target{5}{i_condition}(i_mov) = numel(DAS_all{i_condition}.LT(id_tem1))/numel(DAS_all{i_condition}.LT(id_all));
            end
            test_target{5}{i_condition}(i_mov) = test_target{5}{i_condition}(i_mov)*100;
    end
end
test_p = struct('test_target',[], 'p_all',[]);
test_p.test_target = test_target;
test_p.p_all = zeros(num_condition-1,n_bar);
for i_bar = 1:n_bar
for i_condition = 2:num_condition
   test_p.p_all(i_condition-1,i_bar)=ranksum(test_target{i_bar}{1},test_target{i_bar}{i_condition});
end
end
%================================================================================
save([dir_alt filesep 'DAS_stat_idx.mat'],'DAS_stat', 'idx', 'Cx', 'test_p' ,'-v7.3');
%================================================================================
% else
%     S = load([dir_alt filesep 'DAS_stat_idx.mat']);
%     DAS_stat = S.DAS_stat;
%     idx = S.idx;
end
disp('stat done.')
