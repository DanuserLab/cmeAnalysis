function deleteChannel_Callback(hObject, eventdata, handles)
% Generic callback to be exectuted when a selected channel is removed from
% the graphical settings interface
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

% Get selected properties and returin if empty
selectedProps = get(handles.listbox_selectedChannels, {'String','UserData','Value'});
if isempty(selectedProps{1}) || isempty(selectedProps{3}),return; end

% Delete selected item
% Jenny & Hillary 3/25/24: 
selectedProps{1}(selectedProps{3}) = {''};
selectedProps{2}(selectedProps{3}) = [ ];
set(handles.listbox_selectedChannels, 'String', selectedProps{1});
%check if there is an empty path in selected channels and delete it
%emptyPaths = cellfun("isempty", selectedProps{1});
%set(handles.listbox_selectedChannels, 'String', selectedProps{1}(emptyPaths))
%selectedProps = {selectedProps{1}(~emptyPaths) selectedProps{2} selectedProps{3}};
set(handles.listbox_selectedChannels, 'String', selectedProps{1});
set(handles.listbox_selectedChannels,'UserData',selectedProps{2});
set(handles.listbox_selectedChannels,'Value',max(1,min(selectedProps{3},numel(selectedProps{1}))));
