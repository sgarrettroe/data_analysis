function [out] = fitCLS(dataobj,maxMatrix,flag_plot)
FSlope = fittype( @(m, b, x) m.*x + b, 'coeff', {'m', 'b'}, 'indep', {'x'});

out = struct('fitresult',[],'gof',[],'fitinfo',[]);
for ii = 1:length(dataobj)
    w1 = dataobj(ii).w1;
    y = maxMatrix(ii,:,1)';
    y_std = maxMatrix(ii,:,2)';
    [fitresult,gof,fitinfo] = fit(w1(:), y(:), FSlope, 'StartPoint', [0.5, -20],...
        'lower', [-Inf, -100], 'upper',[10, Inf],'Weight',1./y_std.^2,...
        'Robust','Bisquare');
    
    if flag_plot
        figure(802),clf
        m=15;
        n_contours = 12;   
        ll=linspace(-m,m,n_contours+1);
        map =myMapRGB(n_contours);
        contourf(dataobj(ii).w1,dataobj(ii).w3,dataobj(ii).R)
        colormap(map), hold on,...
        plot(w1, fitresult(w1), 'r','LineWidth',1.5),...
        errorbar(w1,y,y_std,'ro'), hold off;
        xlabel('\omega_1 (cm-1)','FontSize',16)
        ylabel('\omega_3 (cm-1)','FontSize',16)
        set(gcf,'color','white')
        pause
    end
    out(ii).fitresult = fitresult;
    out(ii).gof = gof;
    out(ii).fitinfo = fitinfo;
    out(ii).w1 = w1;
    out(ii).center = y;
    out(ii).center_std = y_std;
end