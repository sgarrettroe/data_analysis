function out = loadReference(fname)
%load reference from the myReference procedure

if exist(fname,'file')==2 %watch out it might be a directory
  temp = load(fname);
elseif exist([fname '.dat'],'file')
  temp = load([fname '.dat']);
  fname = [fname '.dat'];
elseif exist(['Reference_spectrum-' fname],'file')
  temp = load(['Reference_spectrum-' fname]);
  fname = ['Reference_spectrum-' fname];
elseif exist(['Reference_spectrum-' fname '.dat'],'file')
  temp = load(['Reference_spectrum-' fname '.dat']);
  fname = ['Reference_spectrum-' fname '.dat'];
else
  error(['unable to find file ' fname ...
    ' or ' fname '.dat ' ...
    ' or Reference_spectrum-' fname '.dat' ...
    ' in ' pwd]);
end

out = struct('freq',[],...
  'tot',[],...
  'array1',[],...
  'array2',[],...
  'filename',fname);

out.freq = temp(1,2:end);
out.tot = temp(2,2:end);
out.array2 = temp(3,2:end);
out.array1 = temp(2,2:end)+out.array2;
