function s = load2d(basename,pop_time,varargin)
%LOAD2D load 2d data
%
% s = load2d(basename,time) Loads data with a filename that starts with
% <basename> for a population time <time> into a data structure s.
%
% s = load2d(basename,time,'new') Works with the new data format (starting
% some time in 2009???)
%
% s = struct('freq',[],...
%   'time',[],...
%   't2',[],...
%   't3',[],...
%   'w1',[],...
%   'w3',[],...
%   'R',[],...
%   'R1',[],...
%   'R2',[],...
%   'phase',[],...
%   'R1_noise',[],...
%   'R2_noise',[],...  
%   'basename',[],...
%   'undersampling',[],...
%   'centerfreq',[],...
%   'resolution',[],...
%   'zeropad',[],...
%   'time_units',[],...
%   'freq_units',[],...
%   'spec_calib',[],...
%   'pump_probe',[],...
%   'pump_probe_freq',[],...
%   'comment',[]);
%
% call as s = load3d(basename,pop_time,varargin)

s = construct2d;
if nargin == 0
  return
end
s.basename = basename;

flag_new = false;
if ~isempty(varargin)
  if strcmpi(varargin{1},'new')
    flag_new = true;
  end
end

if flag_new
  temp = load([basename,'_T',num2str(pop_time),'.dat']);
else
  temp = load([basename,'_R_T',num2str(pop_time),'.dat']);
end
s.freq = temp(1,2:end);
s.freq_units = 'cm-1';

s.time = temp(2:end,1);
s.time_units = 'fs';

s.R1 = temp(2:end,2:end);

if flag_new
  temp = load([basename,'_T',num2str(pop_time),'_NR.dat']);
else
  temp = load([basename,'_NR_T',num2str(pop_time),'.dat']);
end
s.R2 = temp(2:end,2:end);

s.t2 = pop_time;

%n_delays = length(s.time);
%n_freq = length(s.freq);
%s.R = zeros(n_delays,n_freq);

%s.R1 = zeros(n_delays,n_freq);
%s.R2 = zeros(n_delays,n_freq);
%s.R1_noise = zeros(n_delays,n_delays,n_freq);
%s.R2_noise = zeros(n_delays,n_delays,n_freq);

