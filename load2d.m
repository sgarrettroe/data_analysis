function s = load2d(basename, pop_time, varargin)
% rbLoad2d
% rb, 20110527: wrote function
%
% If there is no varargin given, it will use the naming scheme from before
% summer 2009. If the varargin is 'new', it will use the scheme used
% between summer 2009 and May 2011. If there is a number given, it will use
% the scheme used after May 2011.
%
% input: 
% - basename: the name given in the 2d-mess program
% - pop_time: the population time (T2)
% - time_stamp: the time the measurement started (see note above)
% these three variables are needed to construct the folder and filename
% - scans (opt): default is [0]: import the averaged scan. [n] (where n >
% 0) will import a single specific scan. [n, m, ...] will import multiple
% scans. If noise_flag = true, it will weigh the scans depending on the
% noise
% - noise: will set the noise flag. Default is false. If set to true and
% only 1 scan is imported, it will just import the noise. If multiple scans
% are imported it will weigh the scans based on the noise.
% - meta: will read the metadata. 

s = construct2d;
if nargin == 0
  return
end

% determine type (old, new, newer)
% nothing given: old (before 2009)
% 'new': summer 2009 - may 2011
% 1234: used as timestamp (where 1234 is the time)
type_flag = 0;
if ~isempty(varargin)
  %disp('has varargin')
  if isequal(class(varargin{1}), 'char')  
    %disp('is char')
    if strcmpi(varargin{1},'new')
      %disp('is string new')
      type_flag = 1;
    end
    varargin = varargin(2:end);
  elseif isequal(class(varargin{1}), 'double')
    %disp('is double')
    type_flag = 2;
    time_stamp = varargin{1};
    varargin = varargin(2:end);
  end
end

if type_flag == 2
  % which scans are imported [0]: averaged, [n]: one scan, [n:m] or [n, m, ..]: more scans
  scans = [0];
  % take noise weighted average
  noise_flag = false;
  % import metadata
  meta_flag = true;

  % read arguments
  while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
      case 'scans'
        scans = val;
      case 'noise'
        noise_flag = val;
      case 'meta'
        meta_flag = val;
      otherwise
        error(['load2d: unknown option ',arg])
    end
    varargin = varargin(3:end);
  end


  % fill construct with info
  s.basename = basename;
  s.t2 = pop_time;
  s.time_stamp = time_stamp;

  % construct the base for the name
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
        file_NR_noise = [filebase, '_NR_', num2str(scans(1)), '_noise_.dat'];
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
        file_R_noise = [filebase, '_R_', num2str(i), '_noise.dat'];
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

  % end of importing multiple scans 
  end

  %disp(s)

  if meta_flag 
    file_meta = [folder, filebase, '_meta.txt'];
    s = rbLoadMetadata(s, file_meta);
  end 

  %disp(s)

elseif type_flag == 0 || type_flag == 1
  if type_flag == 1
    temp = load([basename,'_T',num2str(pop_time),'.dat']);
  else
    temp = load([basename,'_R_T',num2str(pop_time),'.dat']);
  end
  
  s.freq = temp(1,2:end);
  s.freq_units = 'cm-1';

  s.time = temp(2:end,1);
  s.time_units = 'fs';
  
  s.R1 = temp(2:end,2:end);

  if type_flag == 1
    temp = load([basename,'_T',num2str(pop_time),'_NR.dat']);
  else
    temp = load([basename,'_NR_T',num2str(pop_time),'.dat']);
  end
  s.R2 = temp(2:end,2:end);

  s.t2 = pop_time;
  
end














