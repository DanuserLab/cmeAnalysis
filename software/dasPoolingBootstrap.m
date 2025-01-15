function [] = dasPoolingBootstrap(dir_DAS, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dir_DAS', @(x) ischar(x));
ip.addParameter('N_bs', 100, @isnumeric);
ip.addParameter('I_sub_samp', 1, @isnumeric);
ip.addParameter('plot_only', false, @islogical);
ip.addParameter('re_write_factor', false, @islogical);
ip.addParameter('N_samp', 20, @isnumeric);
ip.parse(dir_DAS, varargin{:});

dir_DAS = ip.Results.dir_DAS;
N_bs = ip.Results.N_bs;
I_sub_samp = ip.Results.I_sub_samp;
plot_only = ip.Results.plot_only;
re_write_factor = ip.Results.re_write_factor;
N_samp = ip.Results.N_samp;
%======================================================================
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
dir_master = [dir_DAS filesep 'pool'];
dir_data = [dir_master filesep 'data'];
dir_result = [dir_master filesep 'result'];
%======================================================================
S= load([dir_DAS filesep 'DAS_screen.mat']);

con_name_all = [{'siMock'};S.cond_name(2:end)];
date_all = 999999;
for i = 2:max(size(S.cond_date))
    date_all = [date_all, S.cond_date{i}];
end

N=max(size(con_name_all));
n_mov = zeros(N,2);
mat_all = cell(N,1);
for i = 2:N
    mat_all{i} = S.DAS_mat{i,1}{1};
end

clear S;
%======================================================================
%======================================================================
mkdir([dir_data filesep con_name_all{1}]);
mkdir([dir_result filesep con_name_all{1}]);
factor_all = cell(N,1);
A_all = cell(N,1);
for i = 2:N
    factor_all{i} = [1 0];
end
for i = 2:N %1:14
    disp(i)
S_tem = load(mat_all{i});
mkdir([dir_data filesep con_name_all{i}]);
mkdir([dir_result filesep con_name_all{i}]);

S_T = load([S_tem.dir_alt filesep 'ctrl' filesep 'Track_info.mat']);
n_mov(i,1) = max(size(S_T.Track_info));
%==========================================================================
   %-------------------------------------------
   A_all{i} = [];
    for i_mov = 1:n_mov(i,1)
        id_tem = (S_T.Track_info{i_mov}.A>0);
        A_all{i} = [A_all{i};S_T.Track_info{i_mov}.A(id_tem)];
        A_all{i} = randsample(A_all{i},floor(max(size(A_all{i}))*I_sub_samp+0.5));
    end
   %-------------------------------------------
end
save([dir_master filesep 'n_mov.mat'],'n_mov','-v7.3')
disp('----------------data loading done-------------');
%%
if ((exist(fullfile(dir_master, 'factor.mat'), 'file') == 2) == false) || (re_write_factor==true)
A_max_all = 0;
for i = 2:N
    if A_max_all < max(A_all{i})
       A_max_all = max(A_all{i});
    end
end
x = (1:A_max_all);

cdf_emp = cell(N,1);
y = cell(N,1);

cdf_emp{2} = paretotails(A_all{2},0.05,0.95);
y{2} = cdf(cdf_emp{2},x);
for i = 3:N
    cdf_emp{i} = paretotails(A_all{i},0.05,0.95);
    y{i} = cdf(cdf_emp{i},x);
end
disp('----------------initializing for intensity adjustment done-------------');
%%
factor = [ones(N,1), zeros(N,1)];
y_s = y{2};
parfor i = 3:N
    err_min = norm(y{i}-y_s);
    d_factor = [0.1,1];
    max_factor = [20*d_factor(1), 50*d_factor(2)];
    min_factor = [1*d_factor(1), -50*d_factor(2)];
    for factor1 = min_factor(1):d_factor(1):max_factor(1)
        for factor2 = min_factor(2):d_factor(2):max_factor(2)
            cdf_temp = paretotails(A_all{i}*factor1+factor2,0.05,0.95);
            y_temp = cdf(cdf_temp,x);
            err_temp = norm(y_temp-y_s);
            if  err_temp< err_min
                err_min = err_temp;
                factor_all{i} = [factor1,factor2];
            end                
        end
        fprintf('condition# %4.0f is done by %8.3f percent \n',i,factor1/max_factor(1)*100)
    end
end
for i = 3:N
    factor(i,:) = factor_all{i};
end
save([dir_master filesep 'factor.mat'],'factor','-v7.3');
disp('----------------intensity adjustment preparation done-------------');
%%
figure; set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 40, 20]);
  
subplot(1,2,1)
for i = 2:N
    plot(x,y{i},'Linewidth',2)
    hold on
end
set(gca,'Linewidth',2) 
title('before');
xlabel('Intensity')
subplot(1,2,2)
for i = 2:N
    if i == 1
        y_temp = y{i};
    else
    cdf_temp = paretotails(A_all{i}*factor(i,1)+factor(i,2),0.05,0.95);
    y_temp = cdf(cdf_temp,x);
    end
    plot(x,y_temp,'Linewidth',2)
    hold on
end
legend('WT cdf')
set(gca,'Linewidth',2) 
xlabel('Max. Intensity')
title('after');
disp('=====calculation finished=====')
end
%======================================================================
%======================================================================
if plot_only == false
S_tem = load([dir_master filesep 'factor.mat']);
factor = S_tem.factor;
for i = 2:N %1:14
    disp(i)

S_tem = load(mat_all{i});
S_T = load([S_tem.dir_alt filesep 'treated' filesep 'Track_info.mat']);
n_mov(i,2) = max(size(S_T.Track_info));

Track_info = S_T.Track_info;
if i > 2
for i_mov = 1:n_mov(i,2)
    Track_info{i_mov}.A=Track_info{i_mov}.A*factor(i,1)+factor(i,2);
end
end

    mkdir([dir_data filesep 'pooled' filesep 'treated' num2str(i+1)]);
    save([[dir_data filesep 'pooled' filesep 'treated' num2str(i+1)] filesep 'Track_info.mat'],'Track_info','-v7.3');

end
save([dir_master filesep 'n_mov.mat'],'n_mov','-v7.3')
disp('----------------intensity adjustment for treatement conditions done-------------');
%% test2
S_tem = load([dir_master filesep 'n_mov.mat']);
n_mov = S_tem.n_mov;
S_tem = load([dir_master filesep 'factor.mat']);
factor = S_tem.factor;

id_all = []; 
for i = 2:N
id_all = [id_all; [i*ones(n_mov(i,1),1), (1:n_mov(i,1))']];
end
N_id = max(size(id_all));
dir_alt_all = [dir_data filesep 'pooled'];
test_p_all2 = cell(N_bs,1);

for i_bs = 1:N_bs
    disp(i_bs)
Track_info = cell(1,N_samp);  Track_info_alt = cell(1,N_samp);
%id_tem = randsample(N_id,N_samp); id_tem_alt = randsample(N_id,N_samp);
id_tem_all = randsample(N_id,N_samp*2);
id_tem = id_tem_all(1:N_samp); id_tem_alt = id_tem_all(N_samp+1:N_samp*2);
data_all_tem = cell(N+1,1);
for i=2:N
S_tem = load(mat_all{i});
data_all_tem{i+1} = S_tem.data_all{2};
end
data_all_tem{1} = []; data_all_tem{2} = []; 

N_list = id_all(id_tem(:),1);
mov_list = id_all(id_tem(:),2);
N_list_alt = id_all(id_tem_alt(:),1);
mov_list_alt = id_all(id_tem_alt(:),2);
for i = 2:N
    S_tem = load(mat_all{i});
    S_T = load([S_tem.dir_alt filesep 'ctrl' filesep 'Track_info.mat']);
    N_list_tem = 1:N_samp;
    N_list_tem = N_list_tem(N_list==i);
    for j = 1:numel(N_list(N_list==i))
    Track_info{N_list_tem(j)} = S_T.Track_info{mov_list(N_list_tem(j))};
    if N_list_tem(j) > 2
    Track_info{N_list_tem(j)}.A = Track_info{N_list_tem(j)}.A*factor(i,1)+factor(i,2);
    end
    end
    N_list_tem = 1:N_samp;
    N_list_tem = N_list_tem(N_list_alt==i);
    for j = 1:numel(N_list_alt(N_list_alt==i))
    Track_info_alt{N_list_tem(j)} = S_T.Track_info{mov_list_alt(N_list_tem(j))};
    if N_list_tem(j) > 2
    Track_info_alt{N_list_tem(j)}.A = Track_info_alt{N_list_tem(j)}.A*factor(i,1)+factor(i,2);
    end
    end
end

for i_samp = 1:N_samp
    i_N = id_all(id_tem(i_samp),1); i_mov = id_all(id_tem(i_samp),2); 
    i_N_alt = id_all(id_tem_alt(i_samp),1); i_mov_alt = id_all(id_tem_alt(i_samp),2); 
    %------------------------------------------------------
    S_tem = load(mat_all{i_N});
    data_all_tem{1} = [data_all_tem{1}, S_tem.data_all{1}(i_mov)];
    %------------------------------------------------------
    S_tem = load(mat_all{i_N_alt});
    data_all_tem{2} = [data_all_tem{2}, S_tem.data_all{1}(i_mov_alt)];
    %------------------------------------------------------
end

mkdir([dir_data filesep 'pooled' filesep 'ctrl'])
save([[dir_data filesep 'pooled' filesep 'ctrl'] filesep 'Track_info.mat'],'Track_info','-v7.3');
mkdir([dir_data filesep 'pooled' filesep 'treated'])
Track_info=Track_info_alt;
save([[dir_data filesep 'pooled' filesep 'treated'] filesep 'Track_info.mat'],'Track_info','-v7.3');
dir_alt = dir_alt_all;
data_all = data_all_tem;
con_name = [{'Control'}; con_name_all];
save([dir_master filesep 'all.mat'],'data_all','dir_alt','con_name','-v7.3');

    S = load([dir_master filesep 'all.mat']);
    
    pm = dasParameter('data_all',S.data_all,'dir_alt',S.dir_alt,'con_name',S.con_name,...
                  'save_image',false,...
    'plot_fig1', false,'overwriteTrack_info',false,...
    'scheme_p',1,'plot_EpiTIRF',false,'cohort_norm', true,...
    'pdf_conf', false, 'PaperPosition', [0 0 5 5],...
    'plot_visitor',true,'cohort_diff',false,'fig1_mosaic',false,'plot_cohort',false,'fig_disp_mod','none');
dasSingleCondition(pm);
    
    S = load([dir_data filesep 'pooled' filesep 'DAS_stat_idx.mat']);
    test_p_all2{i_bs} = S.test_p;

    save([dir_master filesep 'test_p_all2.mat'],'test_p_all2','-v7.3')
end
disp('========================finished====================')
end
%======================================================================
%======================================================================
S = load([dir_master filesep 'test_p_all2.mat']); 
test_p_all2 = S.test_p_all2;
n_dat = max(size(test_p_all2)); %n_dat=30;
p_all = zeros(N,n_dat,7); perc_d = zeros(N,n_dat,6); fold_d = zeros(N,n_dat,6);  val = zeros(N+1,n_dat,7);
S = load([dir_master filesep 'all.mat']);
con_name = con_name_all;
pm = dasParameter('data_all',S.data_all,'dir_alt',S.dir_alt,'con_name',S.con_name,...
                  'save_image',false,...
    'plot_fig1', false,'overwriteTrack_info',false,...
    'scheme_p',1,'plot_EpiTIRF',false,'cohort_norm', true,...
    'pdf_conf', false, 'PaperPosition', [0 0 5 5],...
    'plot_visitor',true,'cohort_diff',false,'fig1_mosaic',false,'plot_cohort',false,'fig_disp_mod','none');
for i=1:N
    disp(i)
for i_dat=1:n_dat
    for i_bar = 1:6
       p_all(i,i_dat,i_bar) = test_p_all2{i_dat}.p_all(i,i_bar);
       val(i,i_dat,i_bar) = mean(test_p_all2{i_dat}.test_target{i_bar}{i+1});
       %[~,p_all(i,i_dat,i_bar)]=ttest2(test_p_all2{i_dat}.test_target{i_bar}{i+1},test_p_all2{i_dat}.test_target{i_bar}{1});
       perc_d(i,i_dat,i_bar) = 100*(mean(test_p_all2{i_dat}.test_target{i_bar}{i+1})-mean(test_p_all2{i_dat}.test_target{i_bar}{1}))/mean(test_p_all2{i_dat}.test_target{i_bar}{1});
       fold_d(i,i_dat,i_bar) = (mean(test_p_all2{i_dat}.test_target{i_bar}{i+1}))/mean(test_p_all2{i_dat}.test_target{i_bar}{1});
       %perc_d(i,i_dat,i_bar) = mean(test_p_all2{i_dat}.test_target{i_bar}{i+1});
    end
    q_tem1 = test_p_all2{i_dat}.test_target{3}{i+1};
    q_tem2 = test_p_all2{i_dat}.test_target{1}{i+1};
    q_tem3 = test_p_all2{i_dat}.test_target{5}{i+1};
    q_tem1_wt = test_p_all2{i_dat}.test_target{3}{1};
    q_tem2_wt = test_p_all2{i_dat}.test_target{1}{1};
    q_tem3_wt = test_p_all2{i_dat}.test_target{5}{1};
    q_tem  = q_tem1./q_tem2.*q_tem3; q_tem_wt  = q_tem1_wt./q_tem2_wt.*q_tem3_wt;
    p_all(i,i_dat,7) = ranksum(q_tem,q_tem_wt);
    perc_d(i,i_dat,7) = 100*(mean(q_tem)-mean(q_tem_wt))/mean(q_tem_wt);
    %perc_d(i,i_dat,6) = (mean(test_p_all2{i_dat}.test_target{1}{i+1})*mean(test_p_all2{i_dat}.test_target{6}{i+1})-mean(test_p_all2{i_dat}.test_target{1}{1})*mean(test_p_all2{i_dat}.test_target{6}{1}))/(mean(test_p_all2{i_dat}.test_target{1}{1})*mean(test_p_all2{i_dat}.test_target{6}{1}));
    %perc_d(i,i_dat,6) = (1/mean(test_p_all2{i_dat}.test_target{6}{i+1})-1/mean(test_p_all2{i_dat}.test_target{6}{1}))/(1/mean(test_p_all2{i_dat}.test_target{6}{1}));
end
end
%%
i_bar = 1;
i_collect = (1:N); leg_exp = [];
figure;
for i_c = 1:max(size(i_collect))
    %i = i_collect(i_c); histogram(-log10(p_all(i,:,i_bar)),'normalization','pdf','displaystyle','stairs','EdgeColor',pm.col_cond(i,:));
    i = i_collect(i_c); 
    histogram(-log10(p_all(i,:,i_bar)),'normalization','pdf','FaceAlpha',0.5,'FaceColor',pm.col_cond(i,:),'EdgeColor',pm.col_cond(i,:));
    hold on;
    leg_exp = [leg_exp;con_name_all(i)];
end
xlabel('-log_{10}p-value'); ylabel('Prob. density');
h = gca;
plot([-log10(0.001) -log10(0.001)],h.YLim,'--');
leg = legend(leg_exp);legend('boxoff');
leg.ItemTokenSize = leg.ItemTokenSize*0.5;
%%
i_bar_select = [1 4 5];
for i_fig = 1:3
figure; hold on;
i_bar = i_bar_select(i_fig);
FaceAlpha_tem = 0.2;
for i=1:N
    x_to_plot = perc_d(i,:,i_bar); %log2(fold_d(i,:,i_bar))
mean_p = mean(-log10(p_all(i,:,i_bar))); mean_f = mean(x_to_plot); 
notBoxPlot(-log10(p_all(i,:,i_bar)),mean_f);
h=gcf;
h.Children.Children(1).Visible = 'off';
hold on;
for i_ver = 3:4
ver_tem = h.Children.Children(i_ver).Vertices;
cen_tem = mean(ver_tem);
%dx = cen_tem(1) -min(ver_tem(:,1));
dx= 3;
dy = cen_tem(2) -min(ver_tem(:,2));
h.Children.Children(i_ver).Vertices = [cen_tem(1)-dx,cen_tem(2)-dy;...
                                       cen_tem(1)+dx,cen_tem(2)-dy;...
                                       cen_tem(1)+dx,cen_tem(2)+dy;...
                                       cen_tem(1)-dx,cen_tem(2)+dy];
%h.Children.Children(i_ver).FaceAlpha = FaceAlpha_tem;
h.Children.Children(i_ver).FaceColor = pm.col_cond(i,:);
h.Children.Children(i_ver).EdgeColor = [0 0 0];
end
h.Children.Children(2).XData = [cen_tem(1)-dx cen_tem(1)+dx];
h.Children.Children(2).LineWidth = 0.25;
h.Children.Children(3).FaceColor = [0 0 0];
h.Children.Children(4).LineStyle = ':';
h.Children.Children(4).LineWidth = 0.75;

notBoxPlot(x_to_plot,mean_p);
h.Children.Children(2).Color = [1 0 0];
h.Children.Children(2).LineWidth = 0.25;
h.Children.Children(1).Visible = 'off';
x_tem = h.Children.Children(2).XData;
y_tem = h.Children.Children(2).YData;
dx= 0.5;
h.Children.Children(2).XData = [mean(y_tem) mean(y_tem)];
h.Children.Children(2).YData = [mean(x_tem)-dx mean(x_tem)+dx];
for i_ver = 3:4
ver_tem = h.Children.Children(i_ver).Vertices;
cen_tem = mean(ver_tem);
%dx = cen_tem(1) -min(ver_tem(:,1));
dy = cen_tem(2) -min(ver_tem(:,2));
h.Children.Children(i_ver).Vertices = [cen_tem(2)-dy,cen_tem(1)-dx;cen_tem(2)+dy,cen_tem(1)-dx;cen_tem(2)+dy,cen_tem(1)+dx;cen_tem(2)-dy,cen_tem(1)+dx];
h.Children.Children(i_ver).FaceAlpha = FaceAlpha_tem;
h.Children.Children(i_ver).FaceColor = pm.col_cond(i,:);
h.Children.Children(i_ver).EdgeColor = [0 0 0];
end
h.Children.Children(3).FaceAlpha = FaceAlpha_tem;
h.Children.Children(3).FaceColor = [0 0 0];
h.Children.Children(3).FaceAlpha = 1;
h.Children.Children(4).LineStyle = '-';
h.Children.Children(4).EdgeColor = [1 0 1];
end


hold on;
xlim_tem = [-80 40];
plot(xlim_tem, [-log10(0.05) -log10(0.05)],'--'); hold on;
plot(xlim_tem, [-log10(0.01) -log10(0.01)],'--'); hold on;
plot(xlim_tem, [-log10(0.001) -log10(0.001)],'--'); hold on;
xlabel('\Delta_r relative change'); ylabel('-log_{10}p-value'); xlim(xlim_tem);
for i=1:N
    plot(mean(perc_d(i,:,i_bar)),mean(-log10(p_all(i,:,i_bar))),'*','Color',pm.col_cond(i,:));
    hold on;
end
f=get(gca,'Children');
legend(f(N:-1:1),con_name_all);legend('boxoff');
if i_fig == 1
    title('CCP rate (min^{-1}{\mu}m^{-2})')
elseif i_fig == 2
    title('CS init. (min^{-1}{\mu}m^{-2})')
else
    title('CCP%')
end
end


