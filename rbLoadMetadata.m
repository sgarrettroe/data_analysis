function s = rbLoadMetadata(s, filename) %basename, timestamp, pop_time, varargin)


phase_d = 0;
n_shots = 0;
n_scans = 0;

% open the file
fid = fopen(filename, 'r', 'ieee-le.l64', 'windows-932');

tline = fgetl(fid);

if strcmp(tline, 'mess2Dheterodyne_meta_format 1.1')
  while ~feof(fid)
    tline = fgetl(fid);
    
    if strfind(tline, 'Phase') 
      phase = sscanf(tline,'Phase %f?');
      s.phase = (phase)*pi/180;
    end
    
    if strfind(tline, 'Scan')
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      % it prints the start time of the scan, but the last scan is
      % aboprted.
      s.n_scan = str2double(token)-1;
    end

    if strfind(tline, 'Shots') 
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      s.n_shots = str2double(token);
    end

    if strfind(tline, 'T3') 
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      s.t3 = str2double(token);
    end
    
    if strfind(tline, 'Comments') 
      [temp, token] = strtok(tline, ' ');
      %[token, temp] = strtok(token, ' ');
      s.comment = token;
    end
    
  end  
  
  
  
end


fclose(fid);


