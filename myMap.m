function map = myMap(Z,n_levels);
  
  MAX = max(max(Z));
  MIN = min(min(Z));

  %set the value of the deepest black, 
  %   contrast=1 => [0,0,0] 
  %   contrast=.5 =>[0.5 0.5 0.5]
  contrast = .6;
  midtone = contrast/2;
  
  ratio = abs(MAX)/abs(MIN);
  reduce_factor = contrast/(max(MAX,abs(MIN))*2);
    bottom = MIN*reduce_factor;
  top = MAX*reduce_factor;
  range = top-bottom;
  
  step = range/(n_levels);
  offset = 1-midtone;

  v = [bottom:step:top]'+ offset;
%  v=flipud(v);
  map = [v,v,v];
  