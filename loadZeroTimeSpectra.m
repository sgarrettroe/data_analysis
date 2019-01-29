function s = loadZeroTimeSpectra(fname)
%loadZeroTimeSpectra loads the spectra from the zero time determination
%procedure into a structure
if exist(fname,'file')
  fname = strrep(fname,'zt_','spectra_motor');
  dummy1 = load(fname);
  fname = strrep(fname,'spectra_motor','zt_');
  dummy2 = load(fname);
  
else
  error(['Cannot find ',fname,' in path. CWD is ',pwd])
end

freq = dummy1(1,2:end);
time = dummy1(2:end,1)';
spectra = dummy1(2:end,2:end);
zt = dummy2;

s.freq = freq;
s.time = time;
s.spectra = spectra;
s.ztx = zt(:,1);
s.zty = zt(:,2);

figure(1),plot(s.ztx,s.zty,'o')
