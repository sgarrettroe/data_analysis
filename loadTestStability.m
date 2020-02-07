function s = loadTestStability(fname,varargin)

step = 1; %number of beam pairs
t_per_scan = 1;
t_units = 'min';
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case {'step','npairs','n_pairs'}
      step = val;
    case {'time','t','t_per_scan'}
      t_per_scan = val;
    case {'time_units','t_units'}
      t_units = val;
    otherwise
      error(['loadStabilityTest: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

if exist(fname,'file')
  dummy = load(fname);
elseif exist([fname '.dat'],'file')
  dummy = load([fname '.dat']);
else
  error(['Could not find file ' fname ' or ' fname '.dat'])
end
freq = dummy(1,:);
n_freq = length(freq);

scans = cell(1,step);
for i = 1:step
  scans{i} = dummy(1+i:step:end,:);
end

n_scans = size(scans{1},1);
t = (0:n_scans-1)*t_per_scan;

if step >1, warning('the rest of this only works for one pair right now...'),end

figure(1),clf
plot(freq,scans{1}','k')
set(gcf,'windowstyle','docked')

s_hat = fft(scans{1},[],2);
mag = abs(s_hat);
ang = angle(s_hat);
mag(:,1)=0;
mag(:,2)=0;
mag(:,3)=0;
[dummy,ind] = max(mag(1,1:floor(n_freq/2)));
%ind = 5;
dev = (unwrap(ang(:,ind)) - ang(1,5));
%dev = unwrap(dev);
dev =dev*180/pi;
dev3 = dev;
rms_dev3 = std(dev);
figure(4),clf
plot(t,dev3-dev3(1),'ko');
title(['Phase stability of 3D setup ' fname])
ylabel('Phase deviation / deg')
xlabel(['t / ' t_units])
set(gcf,'windowstyle','docked')
