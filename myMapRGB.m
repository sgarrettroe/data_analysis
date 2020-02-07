function map = myMapRGB(n_levels);
%function map = myMapRGB(Z);

%  n_levels = 65;
  n_2 = floor(n_levels/2);
dx = 1/n_2;
  r=[0:dx:1,ones(1,n_2)];
  g=[0:dx:1,1-dx:-dx:0];
  b=[ones(1,n_2),1:-dx:0];
  
  
  map = [r',g',b'];
  
