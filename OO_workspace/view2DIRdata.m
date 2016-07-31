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
N = length(data);
zlimit = 1.0;
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

temp = cropData(data,range1,range3);
f = figure();
clf
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
        %ca
    elseif zlimit > 0 && zlimit <= 1
        [ca, level_list] = myCaxis2(z, n_contours);
        ca = ca * zlimit;
        level_list = level_list * zlimit;
    else
        ca = [-zlimit zlimit];
        level_list = linspace(-zlimit, zlimit, n_contours+2);
    end
    
    caxis(ca)
    if kk == 1
        my2dPlot(x,y,z,'pumpprobe',false,'n_contours',n_contours','zlimit',zlimit);
        ax1 = f.Children(1); % right projection
        ax2 = f.Children(2); % top projection
        ax3 = f.Children(3);
        contourf(ax3,x,y,z,n_contours)
        ax3.Children.LevelList = level_list;
        ax3.XLabel.String = '\omega_1 / 2\pic';
        ax3.YLabel.String = '\omega_3 / 2\pic';
    else
        ax1.Children.XData = sum(z,2);
        ax2.Children.YData = sum(z,1);
        ax3.Children.ZData = z;
        ax3.Children.LevelList = level_list;
    end
        if isfield(data,'scan_number')
            ax2.Title.String = sprintf('t2: %i fs; Run: %03i',data(kk).t2,data(kk).scan_number);
        end
    if kk ~= N
        pause
    end
    %     close(gcf)
end