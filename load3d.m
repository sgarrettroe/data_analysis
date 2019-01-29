function s = load3d(basename)
%LOAD3D load 3d data
% s = struct('freq',[],...
%   'time',[],...
%   't2',[],...
%   't4',[],...
%   't5',[],...
%   'w1',[],...
%   'w3',[],...
%   'w5',[],...
%   'R',[],...
%   'R1',[],...
%   'R2',[],...
%   'R3',[],...
%   'R4',[],...
%   'phase',[],...
%   'R1_noise',[],...
%   'R2_noise',[],...  
%   'R3_noise',[],...  
%   'R4_noise',[],...
%   'basename',basename,...
%   'undersampling',[],...
%   'centerfreq',[],...
%   'resolution',[],...
%   'zeropad',[],...
%   'time_units',[],...
%   'freq_units',[],...
%   'spec_calib',[],...
%   'comment',[]);

s = construct3d;
if nargin == 0
  return
end
s.basename  = basename;

temp = load([basename,'_freq.dat']);
s.freq = temp;
s.freq_units = 'cm-1';

temp = load([basename,'_time.dat']);
s.time = temp;
s.time_units = 'fs';

%if a time file is found
if exist([basename,'_t2t4t5.dat']),
  temp = load([basename,'_t2t4t5.dat']);
  s.t2 = temp(1);
  s.t4 = temp(2);
  s.t5 = temp(3);
end

n_delays = length(s.time);
n_freq = length(s.freq);
s.R = zeros(n_delays,n_delays,n_freq);
s.R1 = zeros(n_delays,n_delays,n_freq);
s.R2 = zeros(n_delays,n_delays,n_freq);
s.R3 = zeros(n_delays,n_delays,n_freq);
s.R4 = zeros(n_delays,n_delays,n_freq);
s.R1_noise = zeros(n_delays,n_delays,n_freq);
s.R2_noise = zeros(n_delays,n_delays,n_freq);
s.R3_noise = zeros(n_delays,n_delays,n_freq);
s.R4_noise = zeros(n_delays,n_delays,n_freq);

temp = load([basename,'.dat']);
%change order to be consistent with Hamm JCP 2006
count = 0;
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R3(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R4(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R1(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R2(j,i,k) = temp(count);
    end
  end
end


temp = load([basename,'_noise.dat']);
count = 0;
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R1_noise(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R2_noise(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R3_noise(j,i,k) = temp(count);
    end
  end
end
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R4_noise(j,i,k) = temp(count);
    end
  end
end


%s.R = s.R1 + s.R2 + s.R3 + s.R4;
