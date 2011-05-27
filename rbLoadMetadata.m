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
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      fprintf(1, 'Phase: %s\n', token);
      %[token idx] = regexp(token, '[0-9]')
      phase_d = str2num(token);
      fprintf(1, 'Ph: %f\n\n', phase_d);
    end
    
    if strfind(tline, 'Scan')
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      %fprintf(1, 'Scan: %s\n', token);
      % it prints the start time of the scan, but the last scan is
      % aboprted.
      s.n_scan = str2double(token)-1;
    end
    
    
    
    if strfind(tline, 'Shots') 
      [temp, token] = strtok(tline, ' ');
      [token, temp] = strtok(token, ' ');
      %fprintf(1, 'Shots: %s\n', token);
      s.n_shots = str2double(token);
    end

    
  end  
  
  %fprintf(1, 'Phase: %.0f\nScans: %.0f\nShots: %.0f\n', phase_d, s.n_scan, s.n_shots)
  
  
end


fclose(fid);


