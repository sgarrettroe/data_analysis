function s = loadSpectrum(fname);
if exist(fname,'file')==2 %careful it could be a directory name
  dummy=load(fname);
elseif exist([fname '.dat'],'file')
  fname = [fname '.dat'];
  dummy=load(fname);
else
  error(['Could not find file ' fname ' or ' fname '.dat ' in pwd]);
end

s = struct('freq',[],...
  'spec',[],...
  'noise',[],...
  'filename',fname);

s.freq = dummy(:,1);
s.spec = dummy(:,2);
if size(dummy,2)>=3
  s.noise = dummy(:,3);
end

