function map = myMapRGB2(n_levels, varargin)
%function map = myMapRGB2(Z);
%paler colors than myMapRGB: min_val = 0.5
%RB, 20110713: now you can set hwo pale you want to have it! Default
%min_val is now 0.2, this gives a bit brighter result.

min_val = 0.2;

while length(varargin) >= 2
  arg = varargin{1};
  val = varargin{2};  
  switch lower(arg)
    case 'brightness'
      min_val = val;
    otherwise
      error(['process2d: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

n_2 = floor(n_levels/2);

dx = (1 - min_val)/n_2;
r=[min_val:dx:1, ones(1,n_2)];
g=[min_val:dx:1, 1-dx:-dx:min_val];
b=[ones(1,n_2), 1:-dx:min_val];

map = [r',g',b'];
  
