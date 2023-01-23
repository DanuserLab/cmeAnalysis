classdef ExternalDetectionProcess < DetectionProcess
    % A concrete class for importing detection results (saved in a .mat file) generated by 3rd party
    % software
    %
    %
    % This class is modified from ExternalSegmentationProcess.m, but
    % specified for detection process.
    %
    % Qiongjing (Jenny) Zou, Feb 2019
%
% Copyright (C) 2023, Danuser Lab - UTSouthwestern 
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
    
    methods(Access = public)
        
        function obj = ExternalDetectionProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            
            super_args{1} = owner;
            super_args{2} = ExternalDetectionProcess.getName();
            super_args{3} = @importExternalDetection;
            if isempty(ip.Results.funParams)
                super_args{4} = ExternalDetectionProcess.getDefaultParams(...
                    owner, ip.Results.outputDir);
            else
                super_args{4} = ip.Results.funParams;
            end
            
            obj = obj@DetectionProcess(super_args{:});
        end

        function sanityCheck(obj)
            sanityCheck@DetectionProcess(obj)

            p = obj.getParameters();
            % Test valid channel index matches input data
            validChannels = find(~cellfun(@isempty, p.InputData));
            assert(isequal(validChannels(:), p.ChannelIndex(:)), 'lccb:set:fatal', ...
                'Selected channels do not match input data\n');

            for i = p.ChannelIndex
                if ~exist(p.InputData{i}, 'dir')
                    error('lccb:set:fatal', ...
                        ['The specified detection directory:\n\n ',p.InputData{i}, ...
                        '\n\ndoes not exist. Please double check your channel path.'])
                end

                if isempty(p.InputData{i})
                    error('lccb:set:fatal', ...
                        ['No proper detection files are detected in:\n\n ',p.InputData{i}, ...
                        '\n\nPlease double check your channel path.'])
                end
            end
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
                 
            if ~obj.owner_.is3D()
                varargout{1} = obj.loadChannelOutput@DetectionProcess(iChan, varargin{:});
            else
                
                % Input check
                outputList = {'movieInfo', 'detect3D','detect3Dall', 'detectionsLabRef'};
                ip = inputParser;
                ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
                ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                    @(x) ismember(x,1:obj.owner_.nFrames_));
                ip.addParameter('useCache',true,@islogical);
                ip.addParameter('iZ',[], @(x) ismember(x,1:obj.owner_.zSize_));
                ip.addParameter('output', outputList{1}, @(x) all(ismember(x,outputList)));
                ip.addParameter('projectionAxis3D','Z', @(x) ismember(x,{'Z','X','Y','three'}));
                ip.parse(iChan, varargin{:})
                output = ip.Results.output;
                iFrame = ip.Results.iFrame;
                projAxis3D = ip.Results.projectionAxis3D;
                iZ = ip.Results.iZ;
                varargout = cell(numel(output), 1);
                ZXRatio = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;
                
                if ischar(output),output={output}; end
                
                for iout = 1:numel(output)
                    switch output{iout}
                        case 'detect3D'
                            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, 'movieInfo');
                            
                            if numel(ip.Results.iFrame)>1
                                v1 = s.movieInfo;
                            else
                                v1 = s.movieInfo(iFrame);
                            end
                            if ~isempty(v1.xCoord) && ~isempty(iZ)
                                % Only show Detections in Z.
                                zThick = 1;
                                tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1), 'VariableNames', {'xCoord','yCoord','zCoord'});
                                valid_states = ((tt.zCoord/ZXRatio)>=(iZ-zThick) & (tt.zCoord/ZXRatio)<=(iZ+zThick));
                                dataOut = tt{valid_states, :};
                                
                                if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                    dataOut = [];
                                end
                            else
                                dataOut = [];
                            end
                            dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                            varargout{iout} = dataOutz;
                            
                        case 'detect3Dall'
                            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, 'movieInfo');
                            
                            if numel(ip.Results.iFrame)>1
                                v1 = s.movieInfo;
                            else
                                v1 = s.movieInfo(iFrame);
                            end
                            if ~isempty(v1.xCoord) && ~isempty(iZ)
                                % Only show Detections in Z.
                                %                             zThick = 1;
                                tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1), 'VariableNames', {'xCoord','yCoord','zCoord'});
                                valid_states = ((tt.zCoord/ZXRatio)>=1 & (tt.zCoord/ZXRatio)<=obj.owner_.zSize_);
                                dataOut = tt{:, :};
                                
                                if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                    dataOut = [];
                                end
                            else
                                dataOut = [];
                            end
                            dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                            varargout{iout} = dataOutz;
                        case 'movieInfo'
                            varargout{iout} = obj.loadChannelOutput@DetectionProcess(iChan, varargin{:});
                        case 'detectionsLabRef'
                            varargout{iout} = load(obj.outFilePaths_{2, iChan}, 'detectionLabRef');
                        otherwise
                            error('Incorrect Output Var type');
                    end
                end
                
            end
            
        end
        
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            if obj.owner_.is3D()
                output(1).name='Detected Objects by zSlice';
                output(1).var = 'detect3D';
                output(1).formatData=@DetectionProcess.formatOutput3D;
                
                output(2) = getDrawableOutput@DetectionProcess(obj);
                output(2).name='Detected Objects';
                output(2).var = 'detect3Dall';
                output(2).formatData=@DetectionProcess.formatOutput3D;
            end
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'External Detection';
        end
        
        function h = GUI()
            h= @externalDetectionProcessGUI;
        end

        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'externalDetection'];
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.InputData = cell(numel(owner.channels_), 1);
        end
    end
end
