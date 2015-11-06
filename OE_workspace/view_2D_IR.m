function view_2D_IR(data,range1,range3)
% Initial data visualization
% This allows us to check and see if the quality of our spectra is good,
% and if our calibration is correct.
temp = cropData(data,range1,range3);
for kk = 1:length(data);
    % display figure
    x = temp(kk).w1;
    y = temp(kk).w3;
    z = temp(kk).R;

    figure(2),clf,
    my2dPlot(x,y,z,'pumpprobe',false,'n_contours',20)
    % add 'zlimit' flag to my2dPlot if you need to zoom in
    pause
end