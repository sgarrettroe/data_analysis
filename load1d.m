function s = load1d(fname)

if ~exist(fname,'file')
  if exist([fname '.dat'])
    fname = [fname '.dat'];
  else
    error(['Could not find file ' fname])
  end
end

%load a 1d spectrum
dummy=load(fname);
s.freq = dummy(:,1)';
s.signal = dummy(:,2)';
s.err = dummy(:,3)';

