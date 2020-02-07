function map = myMapBlueGreen(n_levels);
%function map = myMapRGB(Z);
  
%  n_levels = 65;
  n_2 = (n_levels-1)/2;
dx = 1/n_2;
startat=0.25;
  g=[startat:(1-startat)/n_2:1,ones(1,n_2)];
  r=[startat:(1-startat)/n_2:1,1-dx:-dx:0];
  b=[ones(1,n_2),1:-dx:0];
  
  
  map = [r',g',b'];
  
