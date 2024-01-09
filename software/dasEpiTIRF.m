function [fig_exist,z,z_err] = dasEpiTIRF(pm, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('pm', @(x) isstruct(x));
ip.addParameter('b', 5, @isnumeric);
%==========================================================================
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
ip.parse(pm, varargin{:});
%==========================================================================
pm = ip.Results.pm;
data = pm.data_all;
dir_alt=pm.dir_alt;

b = ip.Results.b;
CohortBounds_s_ccp = pm.CohortBounds_s_ccp;
CohortBounds_s_visitor = pm.CohortBounds_s_visitor;
CohortBounds_s_ftn = pm.CohortBounds_s_ftn;
%==========================================================================
fig_exist = cell(1, 2);
fig_exist{1} = figure;
fig_exist{2} = figure;

for i = 1:2
         figure(fig_exist{i}); 
         set(gcf,'PaperUnits','centimeters');
         set(gcf,'PaperPosition',[0 0 5 5]);
         if i==3
             set(gcf,'PaperPosition',[0 0 11 8]);
         end
         set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
end
%==========================================================================
% new below
%==========================================================================
%==========================================================================
%bound_s = [5 15 25 35 45 55 65 ];
bound_s = CohortBounds_s_ccp;

%bound_s = [0 20 40 60 80 100 120 140 160 180 200];
% fig_tem = dasCohorts(data, DAS_all, idx,dir_alt, pm,...
%     'CohortBounds_s_ccp', CohortBounds_s_ccp,...
%     'CohortBounds_s_visitor', CohortBounds_s_visitor,...
%     'CohortBounds_s_ftn', CohortBounds_s_ftn,...
%     're_read',true,...
%     'MaxIntensityThreshold',0);
% for i_fig = 1:max(size(fig_tem))
%     close(fig_tem{i_fig});
% end
S = load([dir_alt filesep 'DAS_cohort.mat']);

%==========================================================================
n_cond = max(size(data));
z = cell(n_cond,pm.num_clus);
z_err = cell(n_cond,pm.num_clus);
Iepi = cell(n_cond,pm.num_clus);
Itirf = cell(n_cond,pm.num_clus);
n_coh = zeros(n_cond,pm.num_clus);
A_mean_0 = cell(n_cond,pm.num_clus);
A_c_mean_0 = cell(n_cond,pm.num_clus);
for i=1:n_cond
    for ic = 1:pm.num_clus
        n_coh(i,ic) = max(size(S.res{i,ic}.A_mean_t));
        z{i,ic}=cell(n_coh(i,ic));
        z_err{i,ic}=cell(n_coh(i,ic));
        Iepi{i,ic}=cell(n_coh(i,ic));
        Itirf{i,ic}=cell(n_coh(i,ic));
        A_mean_0{i,ic} = zeros(n_coh(i,ic),1);
        A_c_mean_0{i,ic} = zeros(n_coh(i,ic),1);
        for c = 1:n_coh(i,ic)
            A_mean_0{i,ic}(c) = S.res{i,ic}.A_mean_t{c}(2,7);
            A_c_mean_0{i,ic}(c) = S.res_c{i,ic}.A_mean_t{c}(2,7);
        end
    end
end
syms Ie It;
varlist = [Ie It];
n_var = numel(varlist);

for i=1:n_cond
    %figure;hold on;
    for ic = 1:pm.num_clus
        for c = 1:n_coh(i,ic)
            if ((ic == 1) || ((ic == 2) && (CohortBounds_s_visitor(2)>0)) || ((ic == 3) && (CohortBounds_s_ftn(2)>0)))
%-----------------------------------
e=S.res_c{i,ic}.A_mean_t{c}(2,:);
t=S.res{i,ic}.A_mean_t{c}(2,:);
ttem=S.res{i,ic}.A_mean_t{c}(1,:);
%-----------------------------------
t_z0 = b+1; 
ytem=t; 
o_fit = 1;

rsq = zeros(10,1);
for i_fit = 4: 10
t_fit1 = b+1;t_fit2 = b+i_fit;
t_to_fit = [t_fit1:t_fit1+1,t_fit1+3:t_fit2];
[ptem,~]=polyfit(ttem(t_to_fit),ytem(t_to_fit),o_fit);
yfit =  ptem(1) * ttem(t_to_fit) + ptem(2);
yresid = ytem(t_to_fit)-yfit;
SSresid = sum(yresid.^2);
SStotal = (length(ytem(t_to_fit))-1) * var(ytem(t_to_fit));
rsq(i_fit) = 1 - SSresid/SStotal;
end
[~,i_opt] = max(rsq);
t_fit1 = b+1;t_fit2 = b+i_opt;
t_to_fit = [t_fit1:t_fit1+1,t_fit1+3:t_fit2];
[ptem,~]=polyfit(ttem(t_to_fit),ytem(t_to_fit),o_fit);

%plot(ttem(t_fit1:t_fit2),ytem(t_fit1:t_fit2),'*');hold on;plot(ttem(t_fit1):0.1:ttem(t_fit2),polyval(ptem,ttem(t_fit1):0.1:ttem(t_fit2)),'-');
It0 = polyval(ptem,ttem(t_z0));
dptem = polyder(ptem);
kt=polyval(dptem,ttem(t_z0));
%kt=polyval(ptem,ttem(t_fit1+3))-polyval(ptem,ttem(t_fit1+1));
ytem=e;
%t_to_fit = [t_fit1:t_fit2];
rsq = zeros(10,1);
for i_fit = 4: 10
t_fit1 = b+1;t_fit2 = b+i_fit;
t_to_fit = [t_fit1:t_fit1+1,t_fit1+3:t_fit2];
[ptem,~]=polyfit(ttem(t_to_fit),ytem(t_to_fit),o_fit);
yfit =  ptem(1) * ttem(t_to_fit) + ptem(2);
yresid = ytem(t_to_fit)-yfit;
SSresid = sum(yresid.^2);
SStotal = (length(ytem(t_to_fit))-1) * var(ytem(t_to_fit));
rsq(i_fit) = 1 - SSresid/SStotal;
end
[~,i_opt] = max(rsq);
t_fit1 = b+1;t_fit2 = b+i_opt;
t_to_fit = [t_fit1:t_fit1+1,t_fit1+3:t_fit2];
[ptem,~]=polyfit(ttem(t_to_fit),ytem(t_to_fit),o_fit);
%plot(ttem(t_fit1:t_fit2),ytem(t_fit1:t_fit2),'*');hold on;plot(ttem(t_fit1):0.1:ttem(t_fit2),polyval(ptem,ttem(t_fit1):0.1:ttem(t_fit2)),'-');
Ie0 = polyval(ptem,ttem(t_z0));
dptem = polyder(ptem);
ke=polyval(dptem,ttem(t_z0));
%kt=polyval(ptem,ttem(t_fit1+3))-polyval(ptem,ttem(t_fit1+1));
%-----------------------------------
k=kt/ke;
d = -kt/ke*Ie0+It0;
%plot(k*e(6:end-4)+d,'-*');hold on;plot(t(6:end-4),'-*');
DeltaZ = log((k*Ie+d)/It);
sig = vpa(ones(1,n_var));
for i_var = 1:n_var
    sig(i_var) = diff(DeltaZ,varlist(i_var),1);
end
%-----------------------------------
   z{i,ic}{c}=log((k*e+d)./t);

   z_err{i,ic}{c} = zeros(size(z{i,ic}{c}));
   nt = max(size(z{i,ic}{c}));
   [i ic c]
   %-------------------------------------

   for it = 1:nt
       vals = [e(it) t(it)];
       errs = [S.res_c{i,ic}.A_mean_err{c}(it) S.res{i,ic}.A_mean_err{c}(it)];
       err_tem =sqrt((sum((subs(sig,varlist,vals).^2).*(errs.^2))));
       
       z_err{i,ic}{c}(it) = double(err_tem);
   end
            else
                z{i,ic}{c}=[];z_err{i,ic}{c}=[];
            end
        end
    end
end
%==========================================================================
z_max = cell(n_cond,pm.num_clus);
t_max = cell(n_cond,pm.num_clus);
z_mean = cell(n_cond,pm.num_clus);
for i=1:n_cond
    for ic = 1:1:pm.num_clus
        z_max{i,ic} = [];
        z_mean{i,ic} = [];
        for c = 1:n_coh(i,ic)
         [max_tem,t_max_tem] = max(z{i,ic}{c}(b+2:end-b));
         ttem = S.res{i,ic}.A_mean_t{c}(1,b+2:max(size(z{i,ic}{c}))-b);
         z_max{i,ic}=[z_max{i,ic};max_tem];
         t_max{i,ic}=[t_max{i,ic};ttem(t_max_tem)];
         z_mean{i,ic}=[z_mean{i,ic};mean(z{i,ic}{c})];
        end
    end
end
%==========================================================================
for i=1:1
    for ic = 1:2:pm.num_clus
                line_w = 2;line_s='-';
        for c = 1:n_coh(i,ic)                
           z_low = z{i,ic}{c}-z_err{i,ic}{c};
           z_up = z{i,ic}{c}+z_err{i,ic}{c};
           t_tem = S.res{i,ic}.A_mean_t{c}(1,:);
           figure(fig_exist{1}); hold on;
           plot(t_tem(1,b+1:1:end-b),z{i,ic}{c}(b+1:1:end-b),line_s,'linewidth',line_w,'Color',pm.color_clus{ic});
           fill([t_tem(b+1:end-b) t_tem(end-b:-1:b+1)], ...
               [z_low(b+1:end-b) z_up(end-b:-1:b+1)], pm.color_clus{ic},...
                        'EdgeColor', pm.color_clus{ic}, 'FaceAlpha', 0.3,'Parent', gca);
                    if ic == 1
           figure(fig_exist{2}); hold on;
           plot(t_tem(1,b+1:1:end-b),z{i,ic}{c}(b+1:1:end-b),line_s,'linewidth',line_w,'Color',pm.col_cond(i,:));
           fill([t_tem(b+1:end-b) t_tem(end-b:-1:b+1)], ...
               [z_low(b+1:end-b) z_up(end-b:-1:b+1)], pm.col_cond(i,:),...
                        'EdgeColor', pm.col_cond(i,:), 'FaceAlpha', 0.3,'Parent', gca);
                    end
        end
    end
end

for i=1:1
    for ic=1:1
        for c=1:1
%      if ~isempty(t2{i,ic})
%          figure(fig_exist{1});hold on;
%          scatter(t2{i,ic},z{i,ic}{c}(t2{i,ic}+b+1),20,[0 0 0],'filled')
%      end
%      if ~isempty(t3{i,ic})
%          figure(fig_exist{1});hold on;
%          scatter(t3{i,ic},z{i,ic}{c}(t3{i,ic}+b+1),20,[0 0 0],'filled')
%      end
     figure(fig_exist{1});hold on;
it = (S.res{i,ic}.A_mean_t{c}(1,:) == t_max{i,ic}(c));
     scatter(t_max{i,ic}(c),z{i,ic}{c}(it),20,[0 0 0],'filled')
        end
    end
end



for i=2:n_cond
    for ic = 1:1
                line_w = 2;line_s='--';
        for c = 1:n_coh(i,ic)                
           z_low = z{i,ic}{c}-z_err{i,ic}{c};
           z_up = z{i,ic}{c}+z_err{i,ic}{c};
           t_tem = S.res{i,ic}.A_mean_t{c}(1,:);
           figure(fig_exist{2}); hold on;
           plot(t_tem(1,b+1:1:end-b),z{i,ic}{c}(b+1:1:end-b),line_s,'linewidth',line_w,'Color',pm.col_cond(i,:));
           fill([t_tem(b+1:end-b) t_tem(end-b:-1:b+1)], ...
               [z_low(b+1:end-b) z_up(end-b:-1:b+1)], pm.col_cond(i,:),...
                        'EdgeColor', pm.col_cond(i,:), 'FaceAlpha', 0.3,'Parent', gca);
        end
    end
end
figure(fig_exist{1}); hold on;
%ylabel('{\Delta}z/h'); xlabel('t (s)'); 
ylim([-0.2 0.5]); xlim([0 bound_s(end)-b]);
figure(fig_exist{2}); hold on;
%ylabel('{\Delta}z/h'); xlabel('t (s)'); 
ylim([-0.2 0.5]); xlim([0 bound_s(end)-b]);
