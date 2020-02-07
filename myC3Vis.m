function [q135, ah, bh, ch] = myC3Vis(w1,w3,w5,R,varargin)
%myC3Vis.m Plots 3D object and slices to visually gauge the c3 calculation.
%
% call as
%
% [q135, ah, bh, ch]= myC3Vis(w1,w3,w5,R,varargin)
%
% where ah, bh, and ch are the handles to the axes, by row
% and options can be
% 'n_slices'  how many slices in the row
% 'n_contours' how many contours in the contourf plots
% 'spacing' the gap between successive slices
% 'window' fwhm of a windowing function (gaussian)
% 'method' integration method for q135, {|'sum'|,'trapz'}
% 'comment' optional comment string or cell
% 'vis' {|true|,false} whether or not to display the figure
 
%user set parameters
n_slices= 6;
spacing = 35;
n_contours = 10;
flag_boxes = true;
box_color = [0.5 0.5 0.5];
flag_window = false;
window_width = 0;
window = 1; %default, no windowing
method = 'trapz';
comment = [];
flag_vis = true;
flag_scf = true;
flag_auto_mean = true; %caluclate mean freq etc
flag_symmetrize_freq = false;
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case {'n_slices'}
      n_slices = val;
    case {'spacing'}
      spacing = val;
    case {'n_contours'}
      n_contours = val;
    case {'boxes','flag_boxes'}
      flag_boxes = logical(val);
    case 'window'
      flag_window =  true;
      window_width = val;
      %window is calculated below
    case {'method'}
      val = lower(val);
      switch val
        case {'sum','trapz'}
          method = val;
      otherwise
          error('unknown method')
      end
    case 'comment'
      comment = val;
    case 'vis'
      flag_vis = logical(val);
    case 'scf'
       flag_scf = logical(val);
   case 'center_freq'
      flag_auto_mean = false;
       user_mw1 = val(1);
       user_mw3 = val(1);
       user_mw5 = val(1);
    case 'symmetrize_freq'
      flag_symmetrize_freq = val;
    otherwise
      error(['myC3Vis: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

%check if axes are vectors or matrices (built by meshgrid)
%if they are matrices, make them row vectors
if ndims(w1)==3
    w1 = reshape(squeeze(w1(1,:,1)),1,[]);
end
if ndims(w3)==3
    w3 = reshape(squeeze(w3(:,1,1)),1,[]);
end
if ndims(w5)==3
    w5 = reshape(squeeze(w5(1,1,:)),1,[]);
end

%symmetrize freq axes if called for
if flag_symmetrize_freq
  nw1 = length(w1);
  nw3 = length(w3);
  nw5 = length(w5);
  if nw1== nw3 && nw1==nw5
    %if the axes look the same
    ind = find(w1.^2==max(w1.^2));
    if w1(ind)==w3(ind) & w1(ind)==w5(ind)
      % only do this if the value of this point is the same otherwise a
      % more comlpex algorithm would be necessary...
      fprintf('dropping freq points w1 %f w3 %f w5 %f\n',w1(ind),w3(ind),w5(ind));
      w1(ind)=[];
      w3(ind)=[];
      w5(ind)=[];
      R(ind,:,:)=[];
      R(:,ind,:)=[];
      R(:,:,ind)=[];
    end
  else
    error('symmetrize freq not implemented for data with a spectrometer');
  end
end

%Calculate the norm
%originally: norm = sum(R(:));
switch method
  case 'sum'
    norm = sum(R(:));
  case 'trapz'
    %warning('this is not well tested yet...')
    norm = trapz(trapz(trapz(R)));
  otherwise
    error('unknown method');
end
fprintf('norm %f ',norm);

if flag_auto_mean,
  mw1 = sum(w1.*sum(sum(R,3),1))/norm; %sum z then y
  mw3 = sum(w3.*squeeze(sum(sum(R,3),2))')/norm; %sum z then x
  mw5 = sum(w5.*squeeze(sum(sum(R,2),1))')/norm; %sum x then y
  fprintf('mean freq w1 %f w3 %f w5 %f',mw1,mw3,mw5);
else
  mw1 = user_mw1;
  mw3 = user_mw3;
  mw5 = user_mw5;
end
fprintf('\n'); %don't forget newline

%these are now dw1, dw3, dw5
w1 = w1 - mw1;
w3 = w3 - mw3;
w5 = w5 - mw5;

[W1,W3,W5] = meshgrid(w1,w3,w5);
window = ones(size(R));
if flag_window
  window_2sig2 = 2*(window_width/2.355)^2;
  window = exp(-W1.^2/window_2sig2 ...
    -W3.^2/window_2sig2 ...
    -W5.^2/window_2sig2);

  if flag_scf
    %treat the above as a first guess and do at least one iteration
    %correcting the mean frequencies and norm. The norm seems most
    %important...
    old_norm = norm;
    norm = sum(R(:).*window(:));
    mw1 = sum(w1.*sum(sum(R.*window,3),1))/norm; %sum z then y
    mw3 = sum(w3.*squeeze(sum(sum(R.*window,3),2))')/norm; %sum z then x
    mw5 = sum(w5.*squeeze(sum(sum(R.*window,2),1))')/norm; %sum x then y
    
    %       if verbose>=1
    %         disp(['norm ratio ' num2str(norm/old_norm)]);
    %         disp(['mean w1 ' num2str(mw1)]);
    %         disp(['mean w3 ' num2str(mw3)]);
    %         disp(['mean w5 ' num2str(mw5)]);
    %       end
  
    dw1 = w1 - mw1;
    dw3 = w3 - mw3;
    dw5 = w5 - mw5;
      
    [W1,W3,W5]=meshgrid(dw1,dw3,dw5);
  end %end if scf
end %end if window
w3R = W1.*W3.*W5.*R.*window;

%calculate q135
switch method
  case 'sum'
    q135 = sum(sum(sum(W1.*W3.*W5.*R.*window)))/norm;
  case 'trapz'
    %warning('this is not well tested yet...')
    q135 = trapz(trapz(trapz(W1.*W3.*W5.*R.*window)))/norm;
  otherwise
    error('unknown method');
end

if flag_vis==false
  ah = [];
  bh = [];
  ch = [];
  return
end
%disp(['still going and flag_vis = ' num2str(flag_vis)]);

% more or less fixed parameters:
%3d plot params
w=.4;
h=.4;
x1 = 0.5;
y1 = 0.75;

%2d plot params
x2 = 0.5;
x3 = 0.5;
y2 = 0.4;
y3 = 0.15;
%w2d=.2; %default size
h2d=.2; %default size

%start calculating
%for 3d
bx=x1-w/2;
by=y1-h/2;

%size calc
w2dtemp = 0.9/n_slices;
w2d = min([w2dtemp h2d]);
h2d = w2d;

%do the 3d plot
my3dPlot3(w1,w3,w5,R,'n_contours',n_contours)%nb you have to use the same number of contours
ahh = gca;%save axis handle
set(ahh,'Position',[bx by w h])


%level list for w3R
%this ought to be done by myCaxis2... but that touches caxis for the
%3dplot!
 MAX =  max(w3R(:));
 level_list = linspace(-MAX,MAX,n_contours+2);
 dl = level_list(2)-level_list(1); 
 cmin =level_list(1)-dl/2;
 cmax =level_list(end);
%[c2,ll2]=myCaxis2(w3R,n_contours);

w1I = linspace(-spacing*n_slices/2,spacing*n_slices/2,n_slices);
a = zeros(1,n_slices);
b = zeros(1,n_slices);
for i = 1:n_slices,
  bx=x2-w2d*n_slices/2+(i-1)*w2d;
  by=y2-h2d/2;
  
  [W1I,W3I,W5I]=meshgrid(w1I(i),w3,w5);
  sl =interp3(w1,w3,w5,R,W1I,W3I,W5I,'linear',0);%extrap value = 0
  sl = squeeze(sl)';
  a(i)=axes('position',[bx by w2d h2d]);

%  contourf(a(i),w3,w5,sl,n_contours);
  [c,ll]=myCaxis2(sl,n_contours);
  contourf(a(i),w3,w5,sl,ll);
  colormap(myMapRGB2(n_contours));
  caxis(c);%;+(ll(2)-ll(1))/2); %this is still not perfect...
  %box off
  %set(a(i),'XTick',[w1(indw1)],'YTick',[w1(indw1)],'XTickLabel',[],'YTickLabel',[]);
  set(a(i),'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
  %set(a(i),'XTickLabel',[],'YTickLabel',[]);
  line([0 0],[w5(1) w5(end)],'Color',box_color);
  line([w3(1) w3(end)],[0 0],'Color',box_color);
  
  bx=x3-w2d*n_slices/2+(i-1)*w2d;
  by = y3-h2d/2; %0.05; %move the bottom down
  b(i)=axes('position',[bx by w2d h2d]);
  sl =interp3(w1,w3,w5,w3R,W1I,W3I,W5I,'linear',0);%extrap value = 0
  sl = squeeze(sl)';
  contourf(b(i),w3,w5,sl,level_list);
  caxis([cmin cmax]);
%  contourf(b(i),w3,w5,sl,ll2);
%  caxis(c2);
  %box off;
  %set(a(i),'XTick',[w1(indw1)],'YTick',[w1(indw1)],'XTickLabel',[],'YTickLabel',[]);
  set(b(i),'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
  %set(a(i),'XTickLabel',[],'YTickLabel',[]);
  line([0 0],[w5(1) w5(end)],'Color',box_color);
  line([w3(1) w3(end)],[0 0],'Color',box_color);
  
  
end

%do the annotations
arrow_x = [0.05 0.95];
arrow_y = [0.045 0.045];
tic_dy = 0.005;
annotation('doublearrow',arrow_x,arrow_y);
annotation('line',[x3 x3],arrow_y+[-tic_dy, tic_dy]);
text_box_size = 0.05;
text_box_coords = [x3-text_box_size/2, max(arrow_y(1)-tic_dy-text_box_size,0.0001), text_box_size, text_box_size];
string = '0';
annotation('textbox',text_box_coords,...
  'String',string,...
  'HorizontalAlignment','center',...
  'LineStyle','none');

string = num2str(w1I(1));
text_box_coords = [x3-w2d*(n_slices-1)/2-text_box_size/2, max(arrow_y(1)-tic_dy-text_box_size,0.0001), text_box_size, text_box_size];
annotation('textbox',text_box_coords,...
  'String',string,...
  'HorizontalAlignment','center',...
  'LineStyle','none');

string = num2str(w1I(end));
text_box_coords = [x3+w2d*(n_slices-1)/2-text_box_size/2, max(arrow_y(1)-tic_dy-text_box_size,0.0001), text_box_size, text_box_size];
annotation('textbox',text_box_coords,...
  'String',string,...
  'HorizontalAlignment','center',...
  'LineStyle','none');

string = '\omega_1 / 2\pic';
text_box_coords = [0.6, max(arrow_y(1)-tic_dy-text_box_size,0.0001), text_box_size, text_box_size];
annotation('textbox',text_box_coords,...
  'String',string,...
  'HorizontalAlignment','center',...
  'LineStyle','none');

%add comment string here
string = {sprintf('q_{135} = %.2g',q135)};
if flag_window
  string{2} = sprintf('window_{fwhm} = %4.0d cm-1',window_width);
end
if ~isempty(comment)
  if iscell(comment)
    string = [string, comment];
  else
    string = [string, {comment}];
  end
end
text_box_coords = [0.72, 0.8, text_box_size, 0.18];
annotation('textbox',text_box_coords,...
  'String',string,...
  'HorizontalAlignment','left',...
  'LineStyle','none','FitHeightToText','on');

%string = comment
%text_box_coords = [0.6, max(arrow_y(1)-tic_dy-text_box_size,0.0001), text_box_size, text_box_size];
%annotation('textbox',text_box_coords,...
%  'String',string,...
%  'HorizontalAlignment','center',...
%  'LineStyle','none');

set(gcf,'PaperUnits','inches',...
  'PaperPosition',[1 1 6.5 6.5],...
  'PaperOrientation','portrait',...
  'NextPlot','replace',...
  'Position',[1 400 550 550])
figure(gcf)
if nargout>1
  ah = ahh;
  bh = a;
  ch = b;
end
