function view2DIRdata(data,range1,range3,varargin)
% Initial data visualization
% This allows us to check and see if the quality of our spectra is good,
% and if our calibration is correct.
% Current options:
%      'n_contours'
%      'zlimit'

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
for kk = 1:length(data);
    if isempty(temp(kk).R)
        continue
    end
    % display figure
    x = temp(kk).w1;
    y = temp(kk).w3;
    z = temp(kk).R;

    figure(495),clf,
    my2dPlot(x,y,z,'pumpprobe',false,'n_contours',n_contours,'zlimit',zlimit)
    pause
end