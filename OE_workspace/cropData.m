function [dataobj] = cropData(data,range1,range3)

dataobj = struct('w1',[],'w3',[],'R',[]);
for ii = 1:length(data)
    ind1 = find(data(ii).w1>range1(1) & data(ii).w1<range1(2));
    ind3 = find(data(ii).w3>range3(1) & data(ii).w3<range3(2));
    dataobj(ii).w1 = data(ii).w1(ind1);
    dataobj(ii).w3 = data(ii).w3(ind3);
    dataobj(ii).R = data(ii).R(ind3,ind1);
end