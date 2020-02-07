function globalFitPlot(dataMatrix,w1,w3,p_array,gfstruct,varargin)
% default height
height = 135;


while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'displayscans'
            disp_scans = val;
            dataMatrix = dataMatrix(:,:,disp_scans);
            gfstruct.t2_array = gfstruct.t2_array(disp_scans);
            %       case 'yticks'
            %           yticks = val;
            %       case 'xticks'
            %           xticks = val;
        case 'height_px'
            height = val;
        otherwise
            warning(['unknown option ',arg])
    end
    varargin = varargin(3:end);
end
[scale] = gfScale(dataMatrix);
[~,~,n_scans] = size(dataMatrix);
initial = analyticalResponseFunctionsFun(p_array,w1,w3,gfstruct).*scale;
residual = dataMatrix - initial;

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
set(fig,'pos',[x_offset1 y_offset n_scans*w+x_offset1+x_offset2 3*height+1.5*y_offset]);
count = 0;
for ii = 1:n_scans,
    count = count+1;
    
    %data
    had(count) = axes('units','pixels','position',[x_offset1 + (count-1)*w, 2*height + y_offset, w, height]);
    x = w1;
    y = w3;
    z = dataMatrix(:,:,ii);
    [ca, level_list]= myCaxis2(z, n_contours);
    contourf(had(count),x,y,z,level_list);
    axis equal
    line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
    colormap(map)
    caxis(ca);
    t2_lab = text(x(1)+0.1*range(x),y(end)+0.1*(range(y)),[num2str(gfstruct.t2_array(ii)),' ps']); %t2 display
    t2_lab.FontSize = 10; %12;
%     t2_lab.FontWeight = 'bold';
    if count > 1
        set(had(count),'YTick',[])
    end
    set(had(count),'XTickLabel',[])
    
    
    
    
    %analytical response function
    haf(count) = axes('units','pixels','position',[x_offset1 + (count-1)*w, height+y_offset, w, height]);
    z = initial(:,:,ii);
    contourf(haf(count),x,y,z,level_list);
    axis equal
    line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
    colormap(map)
    caxis(ca);
    set(haf(count),'XTickLabel',[])
    if count == 1
        haf(count).YLabel.String = '\omega_{3}/2{\pi}c (cm^{-1})';
        haf(count).YLabel.FontSize = 14;
    end
    if count > 1
        set(haf(count),'YTick',[])
    end
    
    
    
    
    %residual
    har(count) = axes('units','pixels','position',[x_offset1+(count-1)*w, y_offset,w,height]);
    x=w1;
    y=w3;
    z=residual(:,:,ii);
    [ca, level_list]=myCaxis2(z, n_contours);
    contourf(har(count),x,y,z,level_list);
    axis equal
    line([x(1) x(end)],[x(1) x(end)],'Color', [0 0 0]);
    colormap(map)
    caxis(ca);
    if count > 1
        set(har(count),'YTick',[])
    end
    if count == floor(n_scans./2)
        har(count).XLabel.String = '\omega_{1}/2{\pi}c (cm^{-1})';
        har(count).XLabel.FontSize = 14;
    end
end

set(had,'Tickdir','out','ticklength',3.*get(gca,'ticklength'))
set(haf,'Tickdir','out','ticklength',3.*get(gca,'ticklength'))
set(har,'Tickdir','out','ticklength',3.*get(gca,'ticklength'),'XTickLabelRotation',45)