%% Generate centerline slope
% 7 Fitting the Center Line to get Slope
% play with the fitting parameters to get reasonable fitting result
process2dDataInitialize
cd(localdir{1});
process2dDataZhesFcns
% load designated _summary.mat file
file_name = sprintf('%s-%02i_summary.mat',datestring{1},experiment{1});
load(file_name);

clearvars CLS_wt;
slope = zeros(size(data)); % slope for each spectrum
std_slope = zeros(size(data)); %standard deviation of the slope for each spectrum

for ii = 1:length(data);
clearvars out w1 w3
CLS.ind1 = find(data(ii).w1>CLS.range1(1) & data(ii).w1<CLS.range1(2));
CLS.ind3 = find(data(ii).w3>CLS.range3(1) & data(ii).w3<CLS.range3(2));
w1 = data(ii).w1(CLS.ind1);
w3 = data(ii).w3(CLS.ind3);
R = data(ii).R(CLS.ind3,CLS.ind1);
w3a = min(w3):0.1:max(w3);

out = struct('Fit',[]);
% fit each slice of the 2D spectrum (one w1 point, all w3 points) to a
% Gaussian, using the inputs from the external file
    for ij = 1:length(w1)
        out(ij).Fit = fit(w3(:),R(:,ij), F2GausCO2_p, ...
            'StartPoint',CLS.startpoint, 'lower', CLS.lower_bound,...
            'upper',CLS.upper_bound);
        figure(ij + ii - 1),hold on
        plot(w3(:),R(:,ij),'x')
        plot(w3a,out(ij).Fit(w3a));
    end
hold off
% get the peak position and errorbars
std_b = zeros(size(out));
Position = zeros(size(out));
    for m = 1:length(out) 
        dummy = confint(out(m).Fit);
        err_b = dummy(:,2);
        std_b(m) = (err_b(2) - err_b(1))/2; %std of sd
        Position(m) = out(m).Fit.b;
        clear dummy
    end 
    
% start fitting
    FSlope = fittype( @(a, b, x) a.*x + b, 'coeff', {'a', 'b'}, 'indep', {'x'});
    CLS_fit = fit(w1(:), Position(:), FSlope, 'StartPoint', [0.5, -20],...
    'lower', [-Inf, -100], 'upper',[10, Inf],'Weight',1./(std_b).^2);
    data(ii).CLS.fit = CLS_fit;
    data(ii).CLS.w1 = w1;
    data(ii).CLS.Position = Position + d; %d = anharmonicity
    data(ii).CLS.std = std_b;
    slope(ii) = data(ii).CLS.fit.a;
    dummy = confint(data(ii).CLS.fit);
    err_Slope = dummy(:,1);
    std_slope(ii) = (err_Slope(2) - err_Slope(1))/2;

    figure(802),clf
    m=15;
    n_contours = 12;   
    ll=linspace(-m,m,n_contours+1);
    map =myMapRGB(n_contours);
%     contourf(w1,w3,R),...
    contourf(data(ii).w1(ind1),data(ii).w3(ind3),data(ii).R(ind3,ind1))
    colormap(map), hold on,...
    plot(w1, CLS_fit(w1), 'r','LineWidth',1.5),...
    errorbar(w1,Position,std_b,'ro'), hold off;
    xlabel('\omega_1 (cm-1)','FontSize',16)
    ylabel('\omega_3 (cm-1)','FontSize',16)
    set(gcf,'color','white')
%     pause
    clear w1 w3 R
end
CLS.slope = slope;
CLS.std_slope = std_slope;

%%
wt = 1./(CLS.std_slope.^2);
wt_mean = mean(wt);
wt = wt/wt_mean;

CLS.t2 = [data.t2];

clear out
CLS.fit.biexp = fit(CLS.t2(:), CLS.slope(:), F2Exp, 'StartPoint', [1, 1, 2e3, 2e4, 0.05],...
'lower', [0, 0, 0, 0, 0],'Weight',wt);
CLS.fit.monoexp = fit(CLS.t2(:), CLS.slope(:), FExp, 'StartPoint', [1, 10e4, 0],...
'lower', [0, 0, 0],'Weight',wt);
clear wt wt_mean

%%

A = max(CLS.t2);
t = 1:A/100:A;

figure(99),clf
hold on
errorbar(CLS.t2,CLS.slope,CLS.std_slope,'rx');
plot(t,CLS.fit.monoexp(t),'r--');
plot(t,CLS.fit.biexp(t),'b--');
legend('Data','Exponential','Biexponential')
title(data(1).name)
h = gca;
xlabel('t_2 (fs)')
ylabel('CLS (c_2)')
set(h,'Tickdir','out')
hold off
clear A t
%%
beep
flag_save = input('Save output? (1 = "YES"; 0 = "NO"): ');
if flag_save
    save(sprintf('%s-%02i_summary.mat',datestring{1},experiment{1}),...
        'data','CLS','-append')
end
flag_save = 0;