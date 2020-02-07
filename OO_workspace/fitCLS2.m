function [out,c2,c2_std] = fitCLS2(dataobj,maxMatrix,flag_plot)
out = struct('lm',[]);
robustOpts = 'off';

for ii = 1:length(dataobj)
    w1 = dataobj(ii).w1;
    y = maxMatrix(ii,:,1)';
    y_std = maxMatrix(ii,:,2)';
    weights = 1./(y_std.^2);
    
    lm = fitlm(w1,y,'Weights',weights,'RobustOpts',sprintf(robustOpts));
    
    if flag_plot
        % So we can easily plot our line later
        m = table2array(lm.Coefficients(2,1));
        b = table2array(lm.Coefficients(1,1));
        f = @(x) m.*x + b;

        figure(802),clf
        m=15;
        n_contours = 12;   
        ll=linspace(-m,m,n_contours+1);
        map =myMapRGB(n_contours);
        contourf(dataobj(ii).w1,dataobj(ii).w3,dataobj(ii).R)
        colormap(map), hold on,...
        plot(w1,f(w1),'r','linewidth',1.5),...
        errorbar(w1,y,y_std,'ro'), hold off;
        xlabel('\omega_1 (cm-1)','FontSize',16)
        ylabel('\omega_3 (cm-1)','FontSize',16)
        set(gcf,'color','white')
        pause
    end
    out(ii).lm = lm;
    out(ii).w1 = w1;
    out(ii).center = y;
    out(ii).center_std = y_std;
end

c2 = zeros(size(dataobj));
c2_std = c2;
for ii = 1:length(dataobj)
    Z = out(ii).lm.Coefficients(2,1);
    c2(ii) = table2array(Z);
    Z = out(ii).lm.Coefficients(2,2);
    c2_std(ii) = 2*table2array(Z);
end