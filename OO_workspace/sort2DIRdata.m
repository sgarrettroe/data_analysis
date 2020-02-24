function dataOut = sort2DIRdata(dataIn,varargin)
% datastruct = sort2DIRdata(datastruct, varargin) 
% ** WARNING: DOES NOT PROPERLY SORT MULTI-DAY DATA BY RUN NUMBER
% ('scan_number'). **
% SORT2DIRDATA sorts and compresses a 2D-IR data similar to that generated
% by LOAD2DIRDATA. 
%
% using the 'sortby' input will allow you to sort the data by an arbitrary
% field (by default, it sorts by t2 time). e.g.,
%
%      data = sort2DIRdata(data,'sortby','scan_number');
%
% will sort the structure by the field 'scan_number'
%
% Using the 'mode' variable argument will allow you to sort either
% ascending ('ascend') or descending ('descend').
% 
%      data = sort2DIRdata(data,'sortby','scan_number','mode','descend')
%
% Will sort your data by descending order of run numbers ('scan_number')

sort_field = 't2';
sort_mode = 'ascend';
while length(varargin)>=2 %using a named pair
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'mode'
        if strcmp(val,'ascend')||strcmp(val,'descend')
        sort_mode = val;
        else
            warning(['Unknown sorting method: ',val,'. Sorting by ascending order'])
            sort_mode = 'ascend';
        end
    case 'sortby'
        sort_field = val;
    if ~isfield(dataIn,val)
        warning(['Unknown field to sort by',arg,'Sorting by t2'])
        sort_field = 't2';
    end
    otherwise
      warning(['sort2DIRdata: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

% Added a separate function for removing empty elements.
% empty_elems = arrayfun(@(s) all(structfun(@isempty,s)),dataIn);
% temp = dataIn(~empty_elems);
temp = compressStruct(dataIn);
uh = [temp.(sort_field)];
[~,ind] = sort(uh,sort_mode);
dataOut = temp(ind);