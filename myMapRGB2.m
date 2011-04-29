function map = myMapRGB2(n_levels);
%function map = myMapRGB2(Z);
%paler colors than myMapRGB

%  n_levels = 65;
  n_2 = floor(n_levels/2);
  dx = 0.5/n_2;
  r=[0.5:dx:1,ones(1,n_2)];
  g=[0.5:dx:1,1-dx:-dx:0.5];
  b=[ones(1,n_2),1:-dx:0.5];
  
  
  map = [r',g',b'];
  
