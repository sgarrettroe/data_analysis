
for kk = 1:length(data)
    globalfit.ind1 = find(data(kk).w1>globalfit.range1(1) & data(kk).w1<globalfit.range1(2));
    globalfit.ind3 = find(data(kk).w3>globalfit.range3(1) & data(kk).w3<globalfit.range3(2));
    figure(2),clf,
    my2dPlot(data(kk).w1(globalfit.ind1),data(kk).w3(globalfit.ind3),...
        data(kk).R(globalfit.ind3,globalfit.ind1),'pumpprobe',false,...
        'n_contours',20)
    pause
end