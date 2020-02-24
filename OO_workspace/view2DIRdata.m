function view2DIRdata(data,range1,range3,varargin)
% Initial data visualization
% This allows us to check and see if the quality of our spectra is good,
% and if our calibration is correct.
%
% view2DIRdata(data,range1,range3) plots data over ranges of omega1 and
% omega3.
%
% view2DIRdata(data,range1,range3,'name',val)
% Current options:
%      'n_contours'
%      'zlimit'
%
% TODO: Figure out why the color washes out as the amplitude goes down.

N = length(data);
zlimit = 0;
n_contours = 20;

while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'zlimit'
            zlimit = val;
        case 'n_contours'
            n_contours = val;
            if mod(n_contours,2)
                warning('my2dPlot4: Odd number of contour lines may produce unexpected results!')
            end
        otherwise
            warning(['view2DIRdata: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

% Allow the use of empty ranges to plot the whole spectrum
if isempty(range1)
   ind1 = find(~cellfun('isempty',{data.w1}));
   ind1 = ind1(1);
   range1 = [data(ind1).w1(1) data(ind1).w1(end)];
end
if isempty(range3)
   ind3 = find(~cellfun('isempty',{data.w3}));
   ind3 = ind3(1);
   range3 = [data(ind3).w3(1) data(ind3).w3(end)];
end
% 
f = figure();
clf
map=myMapRGB2(n_contours);
colormap(map)
temp = cropData(data,range1,range3);

for kk = 1:N;
    if isempty(temp(kk).R)
        continue
    end
    % display figure
    x = temp(kk).w1;
    y = temp(kk).w3;
    z = temp(kk).R;
    
    if zlimit <= 0
        [ca, level_list]= myCaxis2(z, n_contours);
    elseif zlimit > 0 && zlimit <= 1
        [ca, level_list] = myCaxis2(z, n_contours);
        ca = ca * zlimit;
        level_list = level_list * zlimit;
    else
        ca = [-zlimit zlimit];
        level_list = linspace(-zlimit, zlimit, n_contours+2);
    end
    


    if kk == 1
        my2dPlot(x,y,z,'pumpprobe',false,'n_contours',n_contours','zlimit',zlimit);
        ax1 = f.Children(1); % right projection
        ax2 = f.Children(2); % top projection
        ax3 = f.Children(3);
        ax3.Children(2).LevelList = level_list;
        ax3.XLabel.FontSize = 12;
        ax3.YLabel.FontSize = 12;
    else
        ax1.Children.XData = sum(z,2);
        ax2.Children.YData = sum(z,1);
        ax3.Children(2).ZData = z;
        ax3.Children(2).LevelList = level_list;
        caxis(ax3,ca)
     end
        if isfield(data,'scan_number')
            ax2.Title.String = sprintf('t2: %i fs; Run: %03i',data(kk).t2,data(kk).scan_number);
        end
    if kk ~= N
        pause
    end
end