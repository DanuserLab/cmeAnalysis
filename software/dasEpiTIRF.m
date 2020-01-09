function [fig_exist,z,z_err] = dasEpiTIRF(data, dir_alt, pm, con_name, DAS_stat, DAS_all, idx, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @iscell);
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('con_name', @(x) iscell(x));
ip.addRequired('DAS_stat', @iscell);
ip.addRequired('DAS_all', @iscell);
ip.addRequired('idx', @iscell);
ip.addParameter('b', 5, @isnumeric);
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
ip.parse(data, dir_alt, pm, con_name,DAS_stat, DAS_all, idx, varargin{:});
%==========================================================================
data = ip.Results.data;
pm = ip.Results.pm;
dir_alt=ip.Results.dir_alt;
con_name = ip.Results.con_name;
DAS_stat = ip.Results.DAS_stat;
DAS_all = ip.Results.DAS_all;
idx = ip.Results.idx;
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

t_fit1 = 2+b;t_fit2 = 10+b;
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
    for ic = 1:pm.num_clus
        for c = 1:n_coh(i,ic)
            if ((ic == 1) || ((ic == 2) && (CohortBounds_s_visitor(2)>0)) || ((ic == 3) && (CohortBounds_s_ftn(2)>0)))
%-----------------------------------
e=S.res_c{i,ic}.A_mean_t{c}(2,:);
t=S.res{i,ic}.A_mean_t{c}(2,:);
ttem=S.res{i,ic}.A_mean_t{c}(1,:);
%-----------------------------------
ytem=t;
ptem=polyfit(ttem(t_fit1:t_fit2),ytem(t_fit1:t_fit2),3);
It0 = polyval(ptem,ttem(t_fit1));
dptem = polyder(ptem);
kt=polyval(dptem,ttem(t_fit1));
ytem=e;
ptem=polyfit(ttem(t_fit1:t_fit2),ytem(t_fit1:t_fit2),3);
Ie0 = polyval(ptem,ttem(t_fit1));
dptem = polyder(ptem);
ke=polyval(dptem,ttem(t_fit1));

%-----------------------------------
k=kt/ke;
d = -kt/ke*Ie0+It0;

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
t1 = cell(n_cond,pm.num_clus);
t2 = cell(n_cond,pm.num_clus);
t3 = cell(n_cond,pm.num_clus);

% for i=1:n_cond
%     for ic = 1:4:pm.num_clus
%         z_max{i,ic} = [];
%         z_mean{i,ic} = [];
%         for c = 1:n_coh(i,ic)
%          it = (S.res{i,ic}.A_mean_t{c}(1,:) <= t_max{i,ic}(c)) & (S.res{i,ic}.A_mean_t{c}(1,:) > 1);
%          ztem = z{i,ic}{c}(it);
%          ttem = S.res{i,ic}.A_mean_t{c}(1,it);
%          p = polyfit(ttem,ztem,5);
%          k = polyder(polyder(p));
%          l = (polyder(p));
%          ytem = polyval(p,ttem);
%          ytem1 = polyval(l,ttem);
%          ytem2 = polyval(k,ttem);
%          kappa = abs(ytem2)./(1+ytem1.^2).^1.5;
%          kappa2_tem = [];
%          for item = 2:max(size(ytem2))-1
%              if (kappa(item-1) < kappa(item)) && (kappa(item+1) < kappa(item)) 
%              t2{i,ic} = [t2{i,ic},ttem(item)];
%              kappa2_tem = kappa(item);
%              break;
%              end
%          end
%          
%          if ~isempty(kappa2_tem)
%          for item = 1:max(size(ytem2))-1
%              if (kappa(item) < kappa2_tem) && (kappa(item+1) > kappa2_tem)...
%                      && (ytem2(item) > 0)
%              t3{i,ic} = [t3{i,ic},ttem(item)];
%              break;
%              end
%          end
%          end
%          
%         end
%          if (i==1) && (ic==1)
%          figure(fig_exist{3}); 
%          subplot(1,2,1);plot(ttem,ztem,'o');hold on;plot(ttem,ytem,'-');
%          xlabel('t (s)'); ylabel('{\Delta}z/h'); legend({'{\Delta}z/h','Poly. Fit'},'Location','northwest');   
%          legend('boxoff');
%          subplot(1,2,2);plot(ttem,kappa,'-');
%          xlabel('t (s)'); ylabel('\kappa (s^{-2})');
%          legend({'\kappa'});
%          legend('boxoff');
%          legend('Location','northeast');
%          end
%     end
% end
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
ylim([-0.2 0.7]); xlim([0 bound_s(end)-b]);
figure(fig_exist{2}); hold on;
%ylabel('{\Delta}z/h'); xlabel('t (s)'); 
ylim([-0.2 0.7]); xlim([0 bound_s(end)-b]);
