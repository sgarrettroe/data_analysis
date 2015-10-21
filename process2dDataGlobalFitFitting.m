%%
%tolfun relates to the squared residual in the error
%tolx relates to the displacement in your parameters
tolfun = 1e-17; 
tolx = 1e-10; 
maxfun = 1e3; %maximum number of function evaluations


% set fit parameters
if exist('optimoptions')==2,
    opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
elseif exist('optimset')==2,
%    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
else
    warning('Could not set parameters! Look for the right options functon for your installation.');
end


%We need to set lower (lb) and upper (ub) bounds for our optimization. Be
%careful if you see your fit moving until it sticks to one of these, you
%may need to reassess.

% go!
tic
%fmincon has required parameters of error function and initial guess.
%Documentation has a bunch of additional parameters, most of which we don't
%understand, but the syntax for not using them is to leave them as blanks.
[pfit,err_w] = fmincon(temp.err_fun_w,p0,[],[],[],[],lb,ub,[],opt,globalfit);
toc
globalfit.pfitw = pfit; %save where we finished

for ii = 1:length(pfit)
    fprintf(1,'%20s\t%12f\n',globalfit.pnames{ii},pfit(ii))
end
%%
% calculate the resulting spectra
[temp.fit,temp.extras] = analyticalResponseFunctionsFun(pfit,w1,w3,globalfit);
residual = (temp.data - temp.fit)./temp.sigma;
globalfit.c2 = temp.extras.c2;

% plot result
figure(202),clf
n_scans = length(globalfit.t2_array);
x_offset1 = 0.05;    
x_offset2 = 0.1; %width on the right for parameters
y_offset = 0.05;
w = (1-x_offset1-x_offset2)/n_scans;
height =  0.25;
n_contours = 12;
map = myMapRGB(n_contours);
clear had haf
count = 0;
for ii = 1:length(globalfit.t2_array),
  count = count+1;
  
     %residual
  har(count) = axes('position',[x_offset1 + (count-1)*w, 0.6 + y_offset, w, height]);
  x=w1;
  y=w3;
  
  
  z=temp.residual(:,:,ii);
  [ca, level_list]=myCaxis2(z, n_contours);
  contourf(har(count),x,y,z,level_list);
  line([x(1) x(end)],[x(1) x(end)],'Color', [0 0 0]);
  colormap(map)
  caxis(ca);
  text((x(end)-x(1))*0.1+x(1),(x(end)-x(1))*0.9+x(1),[num2str(globalfit.t2_array(ii)*1000),'fs']); %t2 display
  if count > 1
    set(har(count),'YTick',[])
  end
  
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
  
  z = temp.fit(:,:,ii);
  
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

%parameters
string = [];
for jj = 1:length(globalfit.pnames)
  string = [string,sprintf('%s\t = %7.2f\n',globalfit.pnames{jj},globalfit.pfitw(jj))];
end

h_text =annotation('textbox',[x_offset1+w*(count) y_offset x_offset2 0.5-2*y_offset],'string',string,'Linestyle','none','VerticalAlignment','top');

fprintf(1,'error  \t=%12.3f\n',err_w)
fprintf(1,'error/n\t=%12.3f (ideal = 1)\n',err_w/numel(data));