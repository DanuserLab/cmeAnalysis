function [] = dasMuiltiCondition(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('dir_movie',[],@ischar);
ip.addParameter('dir_DAS',[],@ischar);
ip.addParameter('name_Control','siControl',@ischar);
ip.parse(varargin{:});
%---------------------------------------------------
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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
dir_movie = ip.Results.dir_movie;
dir_DAS = ip.Results.dir_DAS;
dir_rm = {[dir_movie filesep 'siMock']};
name_Control = ip.Results.name_Control;
%=======================================================================
if (exist(fullfile(dir_DAS, 'DAS_screen.mat'), 'file') == 2) == false
M_par=[1.49 108 6.5e-6];
[dir_sub_sub,dir_ctrl_cond,cond_name,cond_date,n_rep_cond,Nf] = getScreenNames(dir_movie,dir_rm,name_Control);
DAS_mat = dir_sub_sub;
for i=2:Nf
    for irep = 1:n_rep_cond(i)
        [DAS_mat_tem] = getDASfiles(dir_DAS,dir_ctrl_cond{i}{irep},dir_sub_sub{i}{irep},cond_name{i},cond_date{i}(irep),M_par);
        DAS_mat{i}{irep} = DAS_mat_tem;
    end
end
clear i; clear irep; clear DAS_mat_tem;
save([dir_DAS filesep 'DAS_screen.mat'],'DAS_mat','cond_name','cond_date','n_rep_cond','Nf');
disp('=============initial preparation done====================');
else
    S = load([dir_DAS filesep 'DAS_screen.mat']);
    DAS_mat = S.DAS_mat; cond_name = S.cond_name; cond_date = S.cond_date; n_rep_cond = S.n_rep_cond;Nf =S.Nf;
end
%=======================================================================
for i=2:Nf
    disp(['loading images for: ' cond_name{i}])
    for irep = 1:n_rep_cond(i)
        S = load(DAS_mat{i}{irep});
        for i_cond = 1:2
        getCellMask(S.data_all{i_cond}, 'Overwrite', false, 'Validate', true);
        close all;
        end
    end
end
disp('=============mask done====================');
%=======================================================================
for i=2:Nf
    disp(['detection/tracking for: ' cond_name{i}])
    for irep = 1:n_rep_cond(i)
        S = load(DAS_mat{i}{irep});
        for i_cond = 1:2
        cmeAnalysis_DTonly(S.data_all{i_cond});
        close all;
        end
    end
end
disp('=============detection&tracking done====================');
%=======================================================================
for i=2:Nf
    for irep = 1:n_rep_cond(i)   
        S = load(DAS_mat{i}{irep});
pdf_conf = false;
save_image = false;
pm = dasParameter('data_all',S.data_all,'dir_alt',S.dir_alt,'con_name',S.con_name,...
                  'save_image',save_image,...
    'plot_fig1', false,'overwriteTrack_info',false,...
    'scheme_p',1,'plot_EpiTIRF',false,'cohort_norm', true,...
    'pdf_conf', pdf_conf, 'PaperPosition', [0 0 5 5],...
    'plot_visitor',true,'cohort_diff',false,'fig1_mosaic',false,'plot_cohort',false,'fig_disp_mod','none');
dasSingleCondition(pm);
clear ans con_name Cx DAS_all DAS_stat data_all dir_alt idx num_condition pm S test_p;
close all;
    end
end
disp('=============DAS done====================');
end
%==========================================================================
%%
function [DAS_mat] = getDASfiles(dir_DAS,dir_ctrl,dir_pert,name_pert,date_pert,M_par)
dir_alt = [dir_DAS filesep name_pert '_' num2str(date_pert)];
DAS_mat = [dir_alt filesep name_pert '_' num2str(date_pert) '.mat'];
DAS_mat_name = [name_pert '_' num2str(date_pert) '.mat'];
%if (exist(fullfile(dir_alt, DAS_mat_name), 'file') == 2) == false
data_all = cell(1,2);
data_all{1} = loadConditionData(dir_ctrl, {''}, {'egfp'},'Parameters',M_par);
data_all{2} = loadConditionData(dir_pert, {''}, {'egfp'},'Parameters',M_par);
for i_cond=1:2
l_tem = [data_all{i_cond}(:).movieLength];
id_tem = l_tem~=mode(l_tem);
data_all{i_cond}(id_tem) = [];
end

con_name = {'siControl',name_pert};
mkdir(dir_alt);
save(DAS_mat,'data_all','con_name','dir_alt');
%end
end
%%
function [dir_sub_sub,dir_ctrl_dcas,dcas_name,dcas_date,n_rep_dcas,Nf] = getScreenNames(dir_MS,dir_rm,dcasCtrl_name)
[dir_sub,dcas_name,~,Nf] = getFolderNames(dir_MS);

N_tem = max(size(dir_rm));
for i_tem = 1:N_tem
for i=1:Nf
   if strcmp(dir_rm{i_tem},dir_sub{i})
       dcas_name(i) = [];
       dir_sub(i) = [];
       Nf = Nf-1;
       break;
   end
end
end
for i=1:Nf
   if strcmp(dcasCtrl_name,dcas_name{i})
       if (i>1) && (i<Nf)
       dcas_name = {dcas_name{i},dcas_name{1:i-1},dcas_name{i+1:end}}';
       dir_sub = {dir_sub{i},dir_sub{1:i-1},dir_sub{i+1:end}}';
       elseif (i==Nf)
       dcas_name = {dcas_name{end},dcas_name{1:end-1}}';
       dir_sub = {dir_sub{end},dir_sub{1:end-1}}';
       end
       break;
   end
end
%==========================================================================
%%
dir_sub_sub = cell(size(dir_sub));
dcas_date = cell(size(dir_sub));
for i=1:Nf
    [dir_tem,~,fdate,~] = getFolderNames(dir_sub{i});
    dir_sub_sub{i}=dir_tem;
    dcas_date{i}=fdate;
end
n_rep_dcas = zeros(Nf,1);
for i=1:Nf
n_rep_dcas(i) = max(size(dir_sub_sub{i}));
end
dir_ctrl_dcas = cell(size(dir_sub_sub));

for i=2:Nf
    for irep = 1:n_rep_dcas(i)
    [dir_ctrl_single] = getCtrlNames(dir_sub_sub{1},dcas_date{1},dcas_date{i}(irep));
    dir_ctrl_dcas{i}(irep) = dir_ctrl_single; 
    end
end
end
%==========================================================================
%%
function [dir_sub,fname,fdate,Nf] = getFolderNames(dir_MS) 
dir_tem = dir(dir_MS);
Ntot = size(dir_tem,1);
Nf = 0;
for i = 1:Ntot
    if dir_tem(i).isdir
        Nf = Nf+1;
    end
end
dir_sub = cell(Nf,1); fname = cell(Nf,1);fdate=zeros(Nf,1);
i_tem=1;
for i = 1:Ntot
    if dir_tem(i).isdir
        fname(i_tem) = {dir_tem(i).name};
        if max(size(dir_tem(i).name)) >= 6
            date_tem = str2num(dir_tem(i).name(1:6));
            if ~isempty(date_tem)
            fdate(i_tem) = date_tem;
            end
        end
        dir_sub(i_tem) = {[dir_MS filesep dir_tem(i).name]}; i_tem = i_tem+1;
    end
end
dir_rm={[dir_MS filesep '.'],...
        [dir_MS filesep '..']};
N_tem = max(size(dir_rm));
for i_tem = 1:N_tem
for i=1:Nf
   if strcmp(dir_rm{i_tem},dir_sub{i})
       fname(i) = [];
       fdate(i) = [];
       dir_sub(i) = [];
       Nf = Nf-1;
       break;
   end
end
end
end
%==========================================================================
%%
function [dir_ctrl_single] = getCtrlNames(dir_ctrl_all,fdate_ctrl,fdate_dcas_single) 
id_tem = fdate_ctrl==fdate_dcas_single;
dir_ctrl_single = dir_ctrl_all(id_tem);
end