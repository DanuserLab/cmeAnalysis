function varargout = listSelectGUI(varargin)
%LISTSELECTGUI launches a GUI to allow selection from a list
%
% SYNOPSIS [selection, selectionList] = listSelectGUI(inputList,maxSelect,moveOrCopy,forceChoice,preselect)
%
% INPUT    inputList   Cell, list of strings or numerical array from which
%                      the user is to select a number of items
%
%          maxSelect   Maximum number of items that can be selected. If
%                      empty, the number of selections is only limited by
%                      the length of the inputList.
%
%          moveOrCopy  optional argument. If 'move', entries are moved from
%                      the left to the right list, if 'copy' (default), the
%                      inputList always remains on the right
%
%          forceChoice optional. If 1 (default), listSelectGUI will launch
%                      on any number of items in the list. If 0 , it will
%                      not open if the number of items is equal to
%                      maxSelect (or 1, if maxSelect is empty)
%
%          preselect   optional. If empty (default), all input is in the
%                      left selection window. A list of indices into
%                      inputList has the selected entries appear in the
%                      right selection window. If indices are negative, it
%                      preselect indicates the list of entries that are not
%                      in the right selection window.
%
% OUTPUT   selection   Indices into inputList of selected items
%
%          selectionList   Cell array containing selected items
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


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @listSelectGUI_OpeningFcn, ...
    'gui_OutputFcn',  @listSelectGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before listSelectGUI is made visible.
function listSelectGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to listSelectGUI (see VARARGIN)

% run opening function until we wait for figure to close.

% if third input argument: remember moveOrCopy (can be 'move' or 'copy')
if length(varargin) > 2 && ~isempty(varargin{3})
    moveOrCopy = varargin{3};
    if ~(strcmpi(moveOrCopy,'move')||strcmpi(moveOrCopy,'copy'))
        % run termination routine
        handles.output = [];
        guidata(hObject,handles)
        h = errordlg('third input argument has to be either ''move'' or '' copy''!');
        uiwait(h)
        close(hObject)
        return
    end
else
    moveOrCopy = 'copy';
end
handles.moveOrCopy = moveOrCopy;

% we allow numbers, strings and cells as input
if isempty(varargin) || isempty(varargin{1})
    % run termination routine
    handles.output = [];
    guidata(hObject,handles)
    h = errordlg('No input list for listSelectGUI');
    uiwait(h)
    close(hObject)
    return
end
inputList = varargin{1};

if iscell(inputList)
    listCell = inputList(:);
elseif ischar(inputList)
    listCell = cellstr(inputList);
elseif isnumeric(inputList)
    listCell = cellstr(num2str(inputList(:)));
else
    % run termination routine
    handles.output = [];
    guidata(hObject,handles)
    h = errordlg('Non-handled input for listSelectGUI');
    uiwait(h)
    close(hObject)
    return
end

% assign list
listLength = length(listCell);
% remember original list
handles.originalList = listCell;

% pre-select
if length(varargin) > 4
    handles.rightIndexList = varargin {5};
    if any(abs(handles.rightIndexList)>listLength)
        error('at least one preselect index is out of range')
    end
    if any(handles.rightIndexList<0)
        handles.rightIndexList = setdiff(1:listLength, abs(handles.rightIndexList));
    end
else
    handles.rightIndexList = [];
end

% store indexLists
if isempty(handles.rightIndexList) || strcmp(moveOrCopy,'copy')
    handles.leftIndexList = (1:listLength)';
else
    handles.leftIndexList = setdiff(1:listLength, handles.rightIndexList);
end
set(handles.LSG_listboxLeft,'String',listCell(handles.leftIndexList));
set(handles.LSG_listboxRight,'String',listCell(handles.rightIndexList));


% the second input argument determines how many items can be selected
if length(varargin) > 1 && ~isempty(varargin{2})
    maxSelect = varargin{2};
    if ~(isnumeric(maxSelect) && isscalar(maxSelect) && all(maxSelect>0))
        % run termination routine
        handles.output = [];
        guidata(hObject,handles)
        h = errordlg('second input argument has to be a scalar >0!');
        uiwait(h)
        close(hObject)
        return
    end
    % set max number of selections
    set(handles.LSG_listboxLeft,'Max',maxSelect);
else
    maxSelect = 0;
    % set max number of selections
    set(handles.LSG_listboxLeft,'Max',listLength);
end

% store maxSelect
handles.maxSelect = maxSelect;


% fourth input argument: forceChoice?
if length(varargin) > 3 && ~isempty(varargin{4})
    forceChoice = varargin{4};
else
    forceChoice = false;
end

% update handles
guidata(hObject,handles);


% check for whether we have to make a choice at all
if ~forceChoice &&listLength <= maxSelect
    % done
else
    % wait for user
    uiwait(handles.listSelectGUI);
end




% --- Outputs from this function are returned to the command line.
function varargout = listSelectGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign varargout if handles still exist
if isempty(handles) %if closed other than by button
    varargout{1} = [];
    varargout{2} = [];
else
    listCell = get(handles.LSG_listboxRight,'String');
    varargout{1} = handles.rightIndexList;
    varargout{2} = listCell;
    set(handles.listSelectGUI,'Visible','off')
    delete(handles.listSelectGUI)
end



% --- Executes on selection change in LSG_listboxLeft.
function LSG_listboxLeft_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_listboxLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if we're getting here through single click: continue.
% on a double click, we'll run moveRight
if strcmp(get(handles.listSelectGUI,'SelectionType'),'normal')
    return
else
    hObject = handles.LSG_moveRight;
    LSG_moveRight_Callback(hObject,[],handles);
end


% --- Executes during object creation, after setting all properties.
function LSG_listboxLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSG_listboxLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LSG_moveRight.
function LSG_moveRight_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_moveRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% we only move right if there is something to be moved
if isempty(get(handles.LSG_listboxLeft,'String'))
    return
end

% get selection, merge with rightIndexList, and then read out from
% originalList
leftSelection = get(handles.LSG_listboxLeft,'Value');
if handles.maxSelect
    % loop through leftSelection. If more than maxSelect items, we throw
    % out the first one
    rightIndexListNew = handles.rightIndexList;
    for i=1:length(leftSelection)
        newItem = handles.leftIndexList(leftSelection(i));
        if any(newItem == rightIndexListNew)
            % we don't add the item
        else
            rightIndexListNew = [rightIndexListNew;...
                newItem];
            if length(rightIndexListNew) > handles.maxSelect
                rightIndexListNew(1) = [];
            end
        end
    end
else
    % just append and sort again
    rightIndexListNew = unique([handles.leftIndexList(leftSelection);...
        handles.rightIndexList]);
end
% show lists
handles.rightIndexList = rightIndexListNew;
set(handles.LSG_listboxRight,'String',...
    handles.originalList(rightIndexListNew));
% set selection to first
set(handles.LSG_listboxLeft,'Value',1);

% if "move", update left list, too
if strcmpi(handles.moveOrCopy,'move')
    handles.leftIndexList(leftSelection) = [];
    if ~isempty(handles.leftIndexList)
        set(handles.LSG_listboxLeft,'String',...
            handles.originalList(handles.leftIndexList));
    else
        % assign empty to listbox for test at top
        set(handles.LSG_listboxLeft,'String','');
    end
end

guidata(hObject,handles)


% --- Executes on button press in LSG_moveLeft.
function LSG_moveLeft_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_moveLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% we only move left if there is something to be moved
if isempty(get(handles.LSG_listboxRight,'String'))
    return
end

% get selection, merge with leftIndexList, and then read out from
% originalList if we're moving. Otherwise, everything is already on the
% left!
rightSelection = get(handles.LSG_listboxRight,'Value');
if strcmpi(handles.moveOrCopy,'move')
    leftIndexListNew = unique([handles.rightIndexList(rightSelection);...
        handles.leftIndexList]);
    handles.leftIndexList = leftIndexListNew;
    set(handles.LSG_listboxLeft,'String',...
        handles.originalList(leftIndexListNew));
end
% set selection to first
set(handles.LSG_listboxRight,'Value',1);

% update right list
handles.rightIndexList(rightSelection) = [];
if ~isempty(handles.rightIndexList)
    set(handles.LSG_listboxRight,'String',...
        handles.originalList(handles.rightIndexList));
else
    % assign empty to listbox for test at top
    set(handles.LSG_listboxRight,'String','');
end

guidata(hObject,handles)

% --- Executes on button press in LSG_Cancel.
function LSG_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% kill figure, which automatically uiresumes
delete(handles.listSelectGUI)

% --- Executes on button press in LSG_OK.
function LSG_OK_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Everything ok -> just uiresume and go into outputFcn
uiresume(handles.listSelectGUI);




% --- Executes on button press in LSG_selAllRight.
function LSG_selAllRight_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_selAllRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set selection to the whole length of the rightIndexList
entriesRight = max(length(handles.rightIndexList),1);
set(handles.LSG_listboxRight,'Value',[1:entriesRight]');


% --- Executes on button press in LSG_selInvLeft.
function LSG_selInvLeft_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_selInvLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% invert selection
newSelection = [1:max(length(handles.leftIndexList),1)]';
newSelection(get(handles.LSG_listboxLeft,'Value')) = [];
if isempty(newSelection)
    newSelection = 1;
end
set(handles.LSG_listboxLeft,'Value',newSelection);


% --- Executes on button press in LSG_selAllLeft.
function LSG_selAllLeft_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_selAllLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set selection to the whole length of the rightIndexList
entriesLeft = max(length(handles.leftIndexList),1);
set(handles.LSG_listboxLeft,'Value',[1:entriesLeft]');


% --- Executes on selection change in LSG_listboxRight.
function LSG_listboxRight_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_listboxRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if we're getting here through single click: continue.
% on a double click, we'll run moveRight
if strcmp(get(handles.listSelectGUI,'SelectionType'),'normal')
    return
else
    hObject = handles.LSG_moveLeft;
    LSG_moveLeft_Callback(hObject,[],handles);
end


% --- Executes during object creation, after setting all properties.
function LSG_listboxRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSG_listboxRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LSG_selInfRight.
function LSG_selInfRight_Callback(hObject, eventdata, handles)
% hObject    handle to LSG_selInfRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% invert selection
newSelection = [1:max(length(handles.rightIndexList),1)]';
newSelection(get(handles.LSG_listboxRight,'Value')) = [];
if isempty(newSelection)
    newSelection = 1;
end
set(handles.LSG_listboxRight,'Value',newSelection);
