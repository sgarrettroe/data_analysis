function s = rbLoad2d(basename, pop_time, time_stamp, varargin)

s = construct2d;
if nargin == 0
  return
end

% which scans are imported [0]: averaged, [n]: one scan, [n..m]: more scans
scans = [0];
% take noise weighted average
noise_flag = false;

while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'scans'
      scans = val;
    case 'noise'
      noise_flag = val;

    otherwise
      error(['rbLoad2d: unknown option ',arg])
  end
  varargin = varargin(3:end);
end



% fill construct with info
s.basename = basename;
s.t2 = pop_time;
s.time_stamp = time_stamp;

folder = [basename, '_', num2str(s.time_stamp), '_T', num2str(s.t2), '/'];
filebase = [basename, '_', num2str(s.time_stamp), '_T', num2str(s.t2)];

% import only a single file
if length(scans) == 1
  
  % import the averaged data
  if scans(1) == 0
    file_R = [filebase, '.dat'];
    file_NR = [filebase, '_NR.dat'];
    if noise_flag == true
      file_R_noise = [filebase, '_R_noise.dat'];
      file_NR_noise = [filebase, '_NR_noise.dat'];      
    end  
  % import a single scan
  else 
    file_R = [filebase, '_R_', num2str(scans(1)), '.dat'];
    file_NR = [filebase, '_NR_', num2str(scans(1)), '.dat'];
    if noise_flag == true
      file_R_noise = [filebase, '_R_', num2str(scans(1)), '_noise.dat'];
      file_NR_noise = [filebase, '_NR_', num2str(scans(1)), '_noise.dat'];
    end
  end
    
  % import R1
  temp = load([folder, file_R]);
  s.R1 = temp(2:end,2:end);

  % extract more info
  s.freq = temp(1,2:end);
  s.freq_units = 'cm-1';
  s.time = temp(2:end,1);
  s.time_units = 'fs';

  % import R2
  temp = load([folder, file_NR]);
  s.R2 = temp(2:end,2:end);
  
  if noise_flag == true
    % import R1 noise
    temp = load([folder, file_R_noise]);
    s.R1_noise = temp(2:end,2:end);  
    % import R2 noise
    temp = load([folder, file_NR_noise]);
    s.R2_noise = temp(2:end,2:end);  
  end
  

% import multiple individual scans
else

  for i = scans
    % construct folder and filename
    file_R = [filebase, '_R_', num2str(i), '.dat'];
    file_NR = [filebase, '_NR_', num2str(i), '.dat'];
    if noise_flag == true
      file_R_noise = [filebase, '_R_', num2str(i), '_noise_.dat'];
      file_NR_noise = [filebase, '_NR_', num2str(i), '_noise_.dat'];
    end
    
    
    temp = load([folder, file_R]);
    
    if noise_flag == false
      % if we shouldn't import noise, set it to ones
      temp_noise = ones(size(temp(:,:)));
    else 
      % otherwise import noise
      temp_noise = load([folder, file_R_noise]);
      % and save it 
      % first check if it exists
      if size(s.R1_noise) == size(temp_noise(2:end,2:end))
        % if so, add it
        s.R1_noise = s.R1_noise + temp_noise(2:end,2:end);
      else 
        % if it doesn't exist, make it
        s.R1_noise = temp_noise(2:end,2:end);
      end  
    end
    
    disp(mean(mean(temp_noise(2:end,2:end))))
    
    % multiply the measurement with the noise (or the ones)
    if size(s.R1) == size(temp(2:end,2:end))
      s.R1 = s.R1 + temp(2:end,2:end) ./ mean(mean(temp_noise(2:end,2:end)));
    else 
      s.R1 = temp(2:end,2:end) ./ mean(mean(temp_noise(2:end,2:end)));
    end    

    if isempty(s.freq) 
      % extract more info
      s.freq = temp(1,2:end);
      s.freq_units = 'cm-1';
      s.time = temp(2:end,1);
      s.time_units = 'fs';
    end
  
    % import R2
    temp = load([folder, file_NR]);
    if noise_flag == false
      % if we shouldn't import noise, set it to ones
      temp_noise = ones(size(temp(:,:)));
    else 
      % otherwise import noise
      temp_noise = load([folder, file_NR_noise]);
      % and save it 
      % first check if it exists
      if size(s.R2_noise) == size(temp_noise(2:end,2:end))
        % if so, add it
        s.R2_noise = s.R2_noise + temp_noise(2:end,2:end);
      else 
        % if it doesn't exist, make it
        s.R2_noise = temp_noise(2:end,2:end);
      end  
    end
    
    disp(mean(mean(temp_noise(2:end,2:end))))
    
    % multiply the measurement with the noise (or the ones)
    if size(s.R2) == size(temp(2:end,2:end))
      s.R2 = s.R2 + temp(2:end,2:end) ./ mean(mean(temp_noise(2:end,2:end)));
    else 
      s.R2 = temp(2:end,2:end) ./ mean(mean(temp_noise(2:end,2:end)));
    end 
  end 
  
  s.R1 = s.R1/length(scans);
  s.R2 = s.R2/length(scans);

  if noise_flag == true
    s.R1_noise = s.R1_noise/length(scans);
    s.R2_noise = s.R2_noise/length(scans);
  end
  
end













