function s = load3dtmp(n_delays,n_freq,basename)
%LOAD3DTMP load 3d temp data
%
%function s = load3dtmp(n_delays,n_freq,basename)
if nargin~=3
  error('wrong number of input arguments call as s = load3dtmp(n_delays,n_freq,basename')
end
s.R = zeros(n_delays,n_delays,n_freq);

temp = load([basename,'.dat']);
count = 0;
for i = 1:n_delays
  for j = 1:n_delays
    for k = 1:n_freq
      count = count+1;
      s.R(j,i,k) = temp(count);
    end
  end
end
