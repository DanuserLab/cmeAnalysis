

function [Dfunc,pdf_t_A, Wnet, w] = dasComputeD(data, Track_info, factor, pm,A_max_all, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x));
ip.addRequired('Track_info', @(x) iscell(x));
ip.addRequired('factor', @(x) isnumeric(x) );
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('A_max_all', @(x) isnumeric(x) );
ip.addParameter('plotD', false, @islogical);
ip.addParameter('data_ratio', 1, @isnumeric);
ip.addParameter('plot_w', false, @islogical);
ip.addParameter('bin_A', 1, @isnumeric);
ip.addParameter('fig_exist',[],@isobject);
ip.addParameter('dplot_lim', [0 150], @isnumeric);
ip.parse(data, Track_info, factor, pm, A_max_all,varargin{:});

bin_A = ip.Results.bin_A;
pm = ip.Results.pm;
factor = ip.Results.factor;
data_ratio = ip.Results.data_ratio;
A_max_all = ip.Results.A_max_all;
tail_perc = 0.05;
if data_ratio < 1
    Track_info = dasRandSamplingSingle(ip.Results.Track_info, data_ratio);
else
    Track_info = ip.Results.Track_info;
end
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
movie_num = size(ip.Results.data,2);
frame_num = ip.Results.data(1).movieLength;
%--------------------------------------------------------------------------
A_bin_num = A_max_all/bin_A;
fig_exist = ip.Results.fig_exist;
dplot_lim = ip.Results.dplot_lim;
%--------------------------------------------------------------------------------------------------------
disp('computing D ...');
z = cell(frame_num,1);
pdf_t_A = zeros(A_bin_num,frame_num);
for i_movie = 1:movie_num
   for i_frame = 1: frame_num-1 
    y = [Track_info{i_movie}.A(Track_info{i_movie}.LT>i_frame,i_frame)*factor(1)+factor(2),Track_info{i_movie}.A(Track_info{i_movie}.LT>i_frame,i_frame+1)*factor(1)+factor(2)];
    z{i_frame} = cat(1,z{i_frame},y);
   end
end

edges = cell(2,1);
edges{1} = 1:bin_A:A_max_all;
edges{2} = 1:bin_A:A_max_all;

%------------------------------------------------
[x1,x2] = meshgrid(edges{1}, edges{2});
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
%------------------------------------------------

w = cell(frame_num,1);
w_norm = cell(frame_num,1);
epsilon = 1;

for i_frame = 1: frame_num-1 

    %fprintf('calculating frame# %d \n',i_frame);

    w{i_frame} = hist3(z{i_frame},'Edges',edges);
    %[f,~] = ksdensity(z{i_frame},xi);
    %f = reshape(f,[A_max_all,A_max_all]);
    %f=f';
    %w{i_frame} = f;
    
    %for iI = 1:A_bin_num
    %test = smoothdata(w{i_frame}(:,iI),'gaussian');
    %w{i_frame}(:,iI) = test;
    %id_temp = w{i_frame}(:,iI) < 0.01*max(w{i_frame}(:,iI));
    %w{i_frame}(id_temp,iI) = 0;
    %end
end




for i_frame = 1: frame_num-1 
   pdf_t_A(:,i_frame) = sum(w{i_frame}(:,:),2);
end

for i_frame = 1: frame_num-1 
    w_norm{i_frame} = w{i_frame};
    for i_A = 1:A_bin_num
      w_norm{i_frame}(i_A,:) = 1 ./ pdf_t_A(i_A,i_frame);
    end
end

Wnet = cell(frame_num,1);
for i_frame = 1: frame_num-1
    wf = w{i_frame};
    wb = w{i_frame}';
    
    %Wnet{i_frame} = log(wf./wb);
    %Wnet{i_frame} = Wnet{i_frame}+log(w_norm{i_frame}./w_norm{i_frame}');
    norm_temp = pdf_t_A(:,i_frame);
    W_temp = wf.*norm_temp'./wb./norm_temp;
    
    %Wnet{i_frame} = dasLogarithm(W_temp);
    Wnet{i_frame} = log(W_temp);
   
   %for iI = 1:A_bin_num
   %test = smoothdata(w{i_frame}(iI,:)/pdf_t_A(iI,i_frame),'gaussian');
   %testT = smoothdata(w{i_frame}(:,iI)'/pdf_t_A(iI,i_frame),'gaussian');
   %Wnet{i_frame}(iI,:) = log(test./testT);
   %end



    
    %Wnet{i_frame} = Wnet{i_frame}./(abs(edges{1}-edges{1}')).^1;
    n_thre = data_ratio*0;
    temp = triu(wf);
    temp = cumsum(temp);
    %temp = (cumsum(temp)./sum(temp));
    %[~,i_start] = min(abs(temp-tail_perc));
    [~,i_start] = min(abs(temp-n_thre));
    temp = triu(wb);
    temp = cumsum(temp);
    %temp = (cumsum(temp)./sum(temp));
    %[~,i_start2] = min(abs(temp-tail_perc));
    [~,i_start2] = min(abs(temp-n_thre));
    for i_A = 1:A_bin_num
        start_temp = max(i_start(i_A),i_start2(i_A));
        norm_temp = numel(pdf_t_A((pdf_t_A(1:start_temp,i_frame)>0),i_frame));
        Wnet{i_frame}(i_A,1:start_temp) = log(sum(wf(i_A,1:start_temp))/sum(wb(1:start_temp,i_A))/pdf_t_A(i_A,i_frame)*mean(pdf_t_A(1:start_temp,i_frame)));
        %Wnet{i_frame}(i_A,1:start_temp-1) = 0;
        %Wnet{i_frame}(i_A,1:start_temp) = log(sum(wf(i_A,1:start_temp))/sum(wb(1:start_temp,i_A))/pdf_t_A(i_A,i_frame)*mean(pdf_t_A((pdf_t_A(1:start_temp,i_frame)>0),i_frame)));
    end
    Wnet{i_frame}(isnan(Wnet{i_frame})) = 0;
    Wnet{i_frame}(isinf(Wnet{i_frame})) = 0;
    
end
%---------------------------
Dfunc = zeros(A_bin_num,frame_num);
for i_frame = 1: frame_num-1 
        wf = w{i_frame};
        wb = w{i_frame}';
    for i_A = 1:A_bin_num
        %i_A_start = i_A-20;
        %if i_A_start <= 0
            i_A_start = 1;
        %end
       %Wnet_temp = Wnet{i_frame}(i_A,i_A_start:i_A);
       %Wnet_temp(Wnet_temp==0) = [];
       %i_A_temp = i_A_start:i_A;
       %i_A_temp_min = min(i_A_temp(Wnet_temp~=0));
       %if ~isempty(i_A_temp_min)
       %    Wnet_temp(1:i_A_temp_min) = [];
       %end
       %norm_temp = numel(Wnet_temp);
       %norm_temp = (i_A-i_A_start+1);
       %norm_temp = 1;
        %Dfunc(i_A,i_frame) = (sum(Wnet{i_frame}(i_A,i_A_start:i_A))  )/(i_A-i_A_start+1);
       num_temp = mean(wf(i_A,1:i_A))/pdf_t_A(i_A,i_frame);
       id_temp = pdf_t_A(1:i_A,i_frame)' >0;
       den_temp = mean(wb(i_A,id_temp)./pdf_t_A(id_temp,i_frame)');
           
       Dfunc(i_A,i_frame) = log( num_temp / den_temp  );

    end
    
    %test = smoothdata(Dfunc(:,i_frame),'gaussian'); 
    %Dfunc(:,i_frame) = test;
        
end
    Dfunc(isnan(Dfunc)) = 0;
    Dfunc(isinf(Dfunc)) = 0;
    
Dfunc_max = max(max(Dfunc));
Dfunc_min = min(min(Dfunc));

%Dfunc(Dfunc>0) = Dfunc(Dfunc>0)/(-Dfunc_min);
%Dfunc = Dfunc/abs(Dfunc_min);

%Dfunc = Dfunc/std(Dfunc(Dfunc~=0))*0.1;

%Dfunc = -Dfunc/mean(Dfunc(Dfunc~=0))*0.02;
%Dfunc = Dfunc/std(std((Dfunc)))*0.1;
%Dfunc(Dfunc > 0) = Dfunc(Dfunc > 0) / mean(Dfunc(Dfunc>0))*0.1;
%Dfunc(Dfunc < 0) = -Dfunc(Dfunc < 0) / mean(Dfunc(Dfunc<0))*0.1;
for i_frame = 1: frame_num-1 
    %test = smoothdata(Dfunc(:,i_frame),'gaussian'); 
 
    %test = smoothdata(Dfunc(:,i_frame),'gaussian'); 
    %Dfunc(:,i_frame) = test;  
end
for i_A = 1:A_bin_num
 
    %test = smoothdata(Dfunc(i_A,:),'gaussian'); 
    %Dfunc(i_A,:) = test;  
end


%Dfunc(Dfunc >0) = 0;
%for i_frame = 1: frame_num-1 
%    Dfunc(:,i_frame) = Dfunc(:,i_frame)*sum(pdf_t_A(:,i_frame))/sum(pdf_t_A(:,1));  
%end

if (ip.Results.plotD == true) 
    if isempty(fig_exist)
        figure();
    else
        figure(fig_exist);
    end
%colormap('bone');
colormap('gray');
imagesc((Dfunc),[-2,2]);
xlim([0 50])
ylim(dplot_lim)
set(gca,'YDir','normal');
colorbar
xlabel('t (s)','interpreter','latex')
ylabel('i (a.u.)','interpreter','latex')
end
%---------------------------
if ip.Results.plot_w == true
figure();
i_plot_w_rep = 15;
colormap('hot')
    imagesc(flipud(Wnet{i_plot_w_rep}),[-10,10]);
    colorbar
    title(num2str(i_plot_w_rep))
    xlim([0 150])
    ylim([A_max_all-150 A_max_all])
end
%---------------------------



