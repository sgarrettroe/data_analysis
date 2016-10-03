function matrixPlot2DIR(dataMatrix,w1,w3,t2_array,dim,varargin)
height = 175;

while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'displayscans'
            disp_scans = val;
            dataMatrix = dataMatrix(:,:,disp_scans);
            t2_array = t2_array(disp_scans);
        case 'height_px'
            height = val;
        otherwise
            warning(['unknown option ',arg])
    end
    varargin = varargin(3:end);
end

n_disp = prod(dim);
m = dim(1);
n = dim(2);

[~,~,n_scans] = size(dataMatrix);

fig = figure();
set(fig,'color','w')
clf
fig.Units = 'pixels';
x_offset1 = 85;
x_offset2 = 50;
y_offset = 75;

w = (range(w1)/range(w3))*height;
n_contours = 14;
map = myMapRGB(n_contours);
set(fig,'pos',[x_offset1 y_offset n*w+x_offset1+x_offset2 m*height+1.5*y_offset]);
count = 0;
for ii = 1:m
    for jj = 1:n
        count = count+1;
        if count > n_scans
            continue
        end
        %data
        had(ii,jj) = axes('units','pixels','position',[x_offset1 + (jj-1)*w, y_offset + (m-1)*height - (ii-1)*height, w, height]);
        x = w1;
        y = w3;
        z = dataMatrix(:,:,count);
        [ca, level_list]= myCaxis2(z, n_contours);
        contourf(had(ii,jj),x,y,z,level_list);
        axis equal
        line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
        colormap(map)
        caxis(ca);
        t2_lab = text(x(1)+0.1*range(x),y(end),[num2str(t2_array(count)),' ps']); %t2 display
        t2_lab.FontSize = 12;
        t2_lab.FontWeight = 'bold';
        
        if ii == ceil(m/2) && jj == 1
            had(ii,jj).YLabel.String = '\omega_{3}/2{\pi}c (cm^{-1})';
            had(ii,jj).YLabel.FontSize = 14;
        end
        if jj ~= 1
            set(had(ii,jj),'YTick',[])
        end
        if count <= n_scans - n
            set(had(ii,jj),'XTickLabel',[]);
        end
        
        if jj == ceil(n./2) && ii == m
            had(ii,jj).XLabel.String = '\omega_{1}/2{\pi}c (cm^{-1})';
            had(ii,jj).XLabel.FontSize = 14;
        end
    end
end
for ii = 1:size(had,1)
    for jj = 1:size(had,2)
        ind(ii,jj) = isa(had(ii,jj),'matlab.graphics.axis.Axes');
    end
end
set(had(ind),'Tickdir','out','ticklength',3.*get(gca,'ticklength'),'XTickLabelRotation',45);