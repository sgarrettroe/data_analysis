
%build initial guess
globalfit.t2_array = t2_array/1000; % fs --> ps
% globalfit.damping = '1exp1fast'; %for Zhe's data
% globalfit.pnames = {'Delta1 (cm-1)','tau1 (ps)',  ...
%     'T2 (ps)','anh (cm-1)','mu12_2','w0 (cm-1)','phi (rad)'};
globalfit.dt = 0.400; %ps
globalfit.n_t = 64; % # of time points
globalfit.w_0_cm = [];% now a fitting parameter -- allows variation of w_0 for the fit
globalfit.bootstrap = []; %[] means no bootstrapping for now

% % Initial guess
% Delta1_cm = 1.5;% wavenumbers
% tau1 = 26; % ps
% T2 = 2.5; % ps
% anh_cm = 24; %anharmonicity
% w0_cm = 2337.5; % wavenumbers
% mu12_2 = 2.2;
% phi = 0;
% 
% p0 = [Delta1_cm tau1 T2 anh_cm mu12_2 w0_cm phi];

disp(p0)

globalfit.p0 = p0; %save where we started

temp.initial = analyticalResponseFunctionsFun(p0,w1,w3,globalfit);
temp.residual = (temp.data - analyticalResponseFunctionsFun(p0,w1,w3,globalfit))./(temp.sigma);

% set up an "error" function
% data is the matrix defined above
temp.err_fun = @(p,globalfit) sum(sum(sum((temp.data - analyticalResponseFunctionsFun(p,w1,w3,globalfit)).^2)));
temp.err_fun_w = @(p,globalfit) sum(sum(sum((temp.data - analyticalResponseFunctionsFun(p,w1,w3,globalfit)).^2./temp.sigma2)));
temp.err_fun_b = @(p,globalfit) sum(sum(sum((temp.data(globalfit.bootstrap) - analyticalResponseFunctionsFun(p,w1,w3,globalfit)).^2./temp.sigma2(globalfit.bootstrap))));

%%
% Plot initial guess and data

figure(201),clf
n_scans = length(globalfit.t2_array);
x_offset1 = 0.05;
y_offset = 0.05;
w = (1-2*x_offset1)/n_scans;
height =  0.25;
n_contours = 14;
map = myMapRGB(n_contours);
clear had haf
count = 0;
for ii = 1:length(globalfit.t2_array),
    count = count+1;
    x=w1;
    y=w3;
    
    %residual
    har(count) = axes('position',[x_offset1 + (count-1)*w, 0.6 + y_offset, w, height]);
    z=temp.residual(:,:,ii);
    [ca,level_list]=myCaxis2(z,n_contours);
    contourf(har(count),x,y,z,level_list);
        line([x(1) x(end)],[x(1) x(end)],'Color', [0 0 0]);
        colormap(map)
        caxis(ca);
        text((x(end)-x(1))*0.1+x(1),(x(end)-x(1))*0.9+x(1),[num2str(globalfit.t2_array(ii)),' ps']); %t2 display
        if count > 1
            set(har(count),'YTick',[])
        end
        set(har(count),'XTickLabel',[])

    %data
    had(count) = axes('position',[x_offset1 + (count-1)*w, 0.3+y_offset, w, height]);
    z = temp.data(:,:,ii);
    [ca, level_list]= myCaxis2(z, n_contours);
    contourf(had(count),x,y,z,level_list);
        line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
        colormap(map)
        caxis(ca);
        text((x(end)-x(1))*0.1+x(1),(x(end)-x(1))*0.9+x(1),[num2str(globalfit.t2_array(ii)*1000),'fs']); %t2 display
        if count > 1
            set(had(count),'YTick',[])
        end
        set(had(count),'XTickLabel',[])

    %fit
    haf(count) = axes('position',[x_offset1 + (count-1)*w, y_offset, w, height]);
    z = temp.initial(:,:,ii);  
    contourf(haf(count),x,y,z,level_list);
        line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
        colormap(map)
        caxis(ca);
        text((x(end)-x(1))*0.1+x(1),(x(end)-x(1))*0.9+x(1),[num2str(globalfit.t2_array(ii)*1000),'fs']); %t2 display
        if count > 1
          set(haf(count),'YTick',[])
        end
end
set(had,'Tickdir','out')
set(haf,'Tickdir','out')