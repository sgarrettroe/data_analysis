function s = average3d(s_in)
%average a bunch of 3d data
%function s_out = average3d(s_in(1:3))
%
% to do: use noise weighted averaging!

s = s_in(1);
n = length(s_in);
R = zeros(size(s_in(1).R1));

for i = 1:n
  if size(s_in(1).R1)~=size(s_in(i).R1), error(['R1s not the same size']),end
  R = R+s_in(i).R1;
end
R = R./n;
s.R1 = R;
for i = 1:n
  if size(s_in(1).R2)~=size(s_in(i).R2), error(['R2s not the same size']),end
  R = R+s_in(i).R2;
end
R = R./n;
s.R2= R;
for i = 1:n
  if size(s_in(1).R3)~=size(s_in(i).R3), error(['R3s not the same size']),end
  R = R+s_in(i).R3;
end
R = R./n;
s.R3 = R;
for i = 1:n
  if size(s_in(1).R4)~=size(s_in(i).R4), error(['R4s not the same size']),end
  R = R+s_in(i).R4;
end
R = R./n;
s.R4 = R;

R = zeros(size(s_in(1).R));
%size(R)
for i = 1:n
  if size(s_in(1).R)~=size(s_in(i).R), error(['Rs not the same size']),end
  %size(s_in(1).R)
  %size(s_in(i).R)
  R = R+s_in(i).R;
end
R = R./n;
s.R = R;

for i = 2:n
  s.basename =[s.basename ' + ' s_in(n).basename];
end
