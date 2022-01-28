% MATLAB-based testing/performance suite for DASC test (part of cmeAnalysis)
% Qiongjing (Jenny) Zou, March 2020
% Test DASC test
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Preconditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_paths = path;
start_dir = pwd;


disp(['PWD:', pwd]);
% Dump path for debugging
s_path = strsplit(path,':');
s_match = (cellfun(@(x) regexp(x,'toolbox'), s_path, 'UniformOutput', false))';
matlab_paths = s_path(cellfun(@isempty, s_match))';
disp('    [MATLAB] current top-level paths....');
disp(matlab_paths);

disp(['Java heap max: ' num2str(java.lang.Runtime.getRuntime.maxMemory/1e9) 'GB'])
disp('Starting DASC test script');


%----Initialization of temp dir

package_name = 'cmeAnalysisDASC';
t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
mkdir(tmpdir);

% ----Gather testing data set
local_test = false;
if local_test
    data_root = 'C:\Users\Jenny\Documents\Data';
    dir_master = [data_root filesep 'DASC_auto_test_lite'];
else
    data_root = '/project/bioinformatics/Danuser_lab/danuser_ci/Data/DASC_auto_test_lite';
    tic;
    disp(['Copying ' data_root 'to ' tmpdir]);
    copyfile(data_root, [tmpdir filesep 'DASC_auto_test_lite']);
    disp('Done copying')
    toc
    
    dir_master = [tmpdir filesep 'DASC_auto_test_lite'];
end

% Analysis Output Directory
saveFolder = fullfile(tmpdir, 'Analysis');
mkdir(saveFolder)
mkdir([saveFolder filesep 'multi_condition/result']);

%% run DASC test script
% adapted from test_wrapper.m by Xinxin Wang with some modifications.

dasMuiltiCondition('dir_movie',[dir_master filesep 'muilti_condition/data'],'dir_DAS',[saveFolder filesep 'multi_condition/result'])
dasPoolingBootstrap([saveFolder filesep 'multi_condition/result'],'I_sub_samp',0.1,'N_samp',5,'N_bs', 1)

%
data_all{1} = loadConditionData([dir_master filesep 'EpiTirf' filesep 'control'], {'TIRF','WT'}, {'egfp','egfp'});
con_name = {'control'};
dir_alt = [saveFolder filesep 'EpiTirf' filesep 'result'];
mkdir(dir_alt)
pm = dasParameter('data_all',data_all,'dir_alt',dir_alt,'con_name',con_name,...
                  'save_image',false,...
    'plot_fig1', false,'overwriteTrack_info',false,...
    'scheme_p',1,'plot_EpiTIRF',true,'cohort_norm', true,...
    'pdf_conf', false, 'PaperPosition', [0 0 5 5],...
    'plot_visitor',true,'cohort_diff',false,'fig1_mosaic',false,'plot_cohort',true,'fig_disp_mod','single');
dasSingleCondition(pm);
close all;
%
data_all = cell(2,1);
data_all{1} = loadConditionData([dir_master '/muilti_condition/data/siControl/190413'], {''}, {'egfp'});
data_all{2} = loadConditionData([dir_master '/muilti_condition/data/siCALM/190413'], {''}, {'egfp'});
con_name = {'siControl','siCALM'};
dir_alt = [saveFolder filesep 'multi_condition' filesep 'result'];
mkdir(dir_alt);
pm = dasParameter('data_all',data_all,'dir_alt',dir_alt,'con_name',con_name,...
                  'save_image',true,...
    'plot_fig1', true,'overwriteTrack_info',false,...
    'scheme_p',1,'plot_EpiTIRF',false,'cohort_norm', true,...
    'pdf_conf', true, 'PaperPosition', [0 0 5 5],...
    'plot_visitor',true,'cohort_diff',false,'fig1_mosaic',false,'plot_cohort',false,'fig_disp_mod','single');
dasSingleCondition(pm);
close all;

disp('Finished DASC test script');

%% Clean up tmpdir
disp('*****************************************');
disp('%% !!!!!!!Cleaning up /tmp/ ');
disp('*****************************************');

cd('~')
ls(tmpdir)
rmdir(tmpdir,'s')
assert(~(exist(tmpdir, 'dir') == 7))

disp('*****************************************');
disp('%% !!!!!!!done cleaning up /tmp/ ');
disp('*****************************************');
