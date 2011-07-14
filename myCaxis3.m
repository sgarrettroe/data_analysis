function varargout = myCaxis2(Z, n_contours)





MAX = max(abs(Z(:)));

level_list = linspace(-MAX, MAX, n_contours);

  level_list = linspace(-MAX,MAX,n_contours+2);
  dl = level_list(2)-level_list(1); % I don't really understand why I have to do this shift
  level_list(1) = level_list(1);%-eps('single');
  level_list(end) = level_list(end);%+eps('single');  
  cmin =level_list(1)-dl/2;
  cmax =level_list(end);%-dl/2+10*eps('single');

  caxis([cmin cmax]);
  
if nargout>=1
  varargout{1} = [cmin cmax];
end
if nargout>=2
  varargout{2} = level_list;
end
