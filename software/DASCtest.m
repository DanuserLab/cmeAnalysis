% MATLAB-based testing/performance suite for DASC test (part of cmeAnalysis)
% Qiongjing (Jenny) Zou, Oct 2019
% Test DASC test
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

package_name = 'cmeAnalysis';
t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
mkdir(tmpdir);

% -------------------------

small_data_set = true

if small_data_set
    disp('Running on truncated data set - small');
    data_root = '/project/bioinformatics/Danuser_lab/danuser_ci/Data/DASC_test/mini_test'
    tic;
    disp(['Copying ' data_root 'to ' tmpdir]);
    copyfile(data_root, [tmpdir filesep 'mini_test']);
    disp('Done copying')
    toc
    matPath = fullfile(tmpdir, 'mini_test', 'test_input.mat');
else
    disp('Running on large data set - may take some time...');
    data_root = '/project/bioinformatics/Danuser_lab/danuser_ci/Data/DASC_test'
    disp('Will not copy data to tmpdir, since data is too big...');
    matPath = fullfile(data_root, 'test_input.mat');
end

% data_root = '/project/bioinformatics/Danuser_lab/danuser_ci/Data/DASC_test_wTiffs'
% data_path1 = fullfile(data_root, 'control');
% data_all{1} = loadConditionData(data_path1, {''}, {'eGFP'});
% data_path2 = fullfile(data_root, 'siCALM');
% data_all{2} = loadConditionData(data_path2, {''}, {'eGFP'});
load(matPath);
dir_alt = fullfile(tmpdir, 'Analysis');
mkdir(dir_alt)
con_name = {'siControl','siCALM'};
pm = dasParameter(size(con_name,2),dir_alt);

%% run DASC test script
tic
wrapper() % this a test script for this DASC test
toc
close all
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




