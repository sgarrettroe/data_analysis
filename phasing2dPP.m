function [phase,t0_bin_shift,analysis]=phasing2dPP(time,igram)
%phasing2dPP Do the phasing data analysis for the pump-probe 2d
%spectrometer


t = time;
n_pairs = size(igram,1);

norm = zeros(1,n_pairs);
for i = 1:n_pairs
  norm(i) = max(igram(i,:));
  igram(i,:) = igram(i,:)./norm(i);
end

%
% process scans for fft etc
%
[abs_array,phase_array,w,i_max] = phasingProcessScans(t,igram);

%
% calculate the final phase and the index of the time closest to 0
%
[phase,ph,t0_bin_shift,i_fit]  = phasingFinalPhase(n_pairs,abs_array,phase_array,w,i_max);

%
% plot interferograms and ffts
%
phasingPlot(n_pairs,t,igram,w,abs_array,phase_array,i_max,i_fit);

%
%    output results
%

 analysis.w = w;
 analysis.abs = abs_array;
 analysis.phase = phase_array;
 analysis.ph = ph;
 analysis.i_max = i_max;
 analysis.t0_bin_shift = t0_bin_shift;

%----------------------------------------------------------------------
% 
% subfunctions
%
function [abs_array,phase_array,w,i_max] = phasingProcessScans(t,igram);
n_t = length(t);

%index to the first guess for t0
[~,i_t0_initial] = min(t.^2);

%shift time points that are less than 0 to the left circularly
igram = circshift(igram,[0 -(i_t0_initial-1)]);

igram = igram - mean(igram,2);

%calculate the fft
%abs_array = fft(fftshift(igram,2)')';
abs_array = fft(igram);
phase_array = angle(abs_array);
abs_array = abs(abs_array);

%generate the freq axis
w = fftFreqAxis(t,'time_units','fs','shift','off');

%take only half the data
ind = floor(n_t/2);
abs_array = abs_array(:,1:ind);
phase_array = phase_array(:,1:ind);
w = w(1:ind);
%find max of the spectrum from pair 1/2 (strong)
[~,i_max]=max(abs_array(1,:));

%---------------------------------------------------------------------
%
% calculate the final phase
%
function [phase,ph,delta_t_fringes,i_fit] = phasingFinalPhase(n_pairs,abs_array,phase_array,w,i_max)
%phasingFinalPhase
global c_cmfs wavenumbersToInvFs fringeToFs
n_fit_points = 5;
fit_points_offset = (n_fit_points-1)/2;

%rough guess
%tau = 1/w(i_max)/wavenumbersToInvFs;
disp(['rough guess frequency ' num2str(w(i_max))]);

w0 = peakpos(w(i_max-2:i_max+2),abs_array(1,i_max-2:i_max+2));
disp(['refined guess frequency ' num2str(w0)]);
%refine guess
tau = 1/w0/wavenumbersToInvFs;

i_fit = i_max-fit_points_offset:i_max+fit_points_offset; %take n_points points around the max

%this is the phase of the interferogram (mean value)
ph = zeros(1,n_pairs);
switch n_pairs
  case 1
    ph(1) = mean(unwrap(phase_array(1,i_fit)))*180/pi;
    simple_phase = ph*pi/180;
  case 2
    ph(1) = mean(unwrap(phase_array(1,i_fit)))*180/pi;
    ph(2) = mean(unwrap(phase_array(2,i_fit)))*180/pi;
    simple_phase = mean(unwrap(phase_array(1,i_fit))) ...
      -mean(unwrap(phase_array(2,i_fit)));
  case 3
    ph(1) = mean(unwrap(phase_array(1,i_fit)))*180/pi;
    ph(2) = mean(unwrap(phase_array(2,i_fit)))*180/pi;
    ph(3) = mean(unwrap(phase_array(3,i_fit)))*180/pi;
    
    simple_phase = mean(unwrap(phase_array(1,i_fit)))...
      +mean(unwrap(phase_array(2,i_fit))) ...
      -mean(unwrap(phase_array(3,i_fit)));
end
disp(['without moving motors the phase should be ',...
  num2str(rem(simple_phase*180/pi,360)),''''])

%this is the shift of the pulse envelope calculated from the slope of the
%phase
delta_t_fs = zeros(1,n_pairs);
delta_t_fringes = zeros(1,n_pairs);
shift_phase= zeros(1,n_pairs);
for i = 1:n_pairs
  %phase_array(i,i_fit)
  %unwrap(phase_array(i,i_fit))
  p = polyfit(w(i_fit),unwrap(phase_array(i,i_fit)),1);
  dph_dnu = p(1);
  delta_t_fs(i) = dph_dnu/(2*pi*c_cmfs);
  delta_t_fringes(i) = round(delta_t_fs(i)/fringeToFs);
  shift_phase(i) = delta_t_fringes(i)*fringeToFs/tau*360;
end
string = sprintf('pair\tdt_fs\t\tdt_fringes');
disp(string)
for i = 1:n_pairs
  string = sprintf('%i\t\t%6.1f\t\t%4i\t%6.1f',i,delta_t_fs(i),delta_t_fringes(i),shift_phase(i));
  disp(string)
end

%this is my attempt to bring those numbers together to make the final phase
switch n_pairs
  case 1 %2d pp case
    phase = ph(1) - shift_phase(1);
  case 2 %2d case
%     disp(['pair 1 ' num2str(ph1)])
%     disp(['pair 2 ' num2str(ph2)])
%     disp(['diff  ' num2str(ph1-ph2)])
%     disp(['2*ph2 ' num2str(ph1-2*ph2)])

    %rem_phase = shift_phase(1)-shift_phase(2);
    %phase = simple_phase*180/pi - rem_phase;
    % or 
    phase = ph(1)-shift_phase(1)-ph(2)+shift_phase(2);
    phase_alt = ph(1)-shift_phase(1) - 2*(ph(2)-shift_phase(2));
    phase_alt2 = ph(1)-shift_phase(1) - 2*ph(2)+shift_phase(2);
  case 3
    rem_phase = shift_phase(1)+shift_phase(2)-shift_phase(3);
    phase = simple_phase*180/pi - rem_phase;    
  otherwise
    error(['n_pairs must be 1, 2 or 3! it is n_pairs = ',num2str(n_pairs)]);
end
% 
phase = rem(phase,360);
string = sprintf('Final phase is %6.1f',phase);
disp('****************************')

%debugging phasing
switch n_pairs
  case 2
    disp(string)
    disp('****************************')
    string = sprintf('%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f',...
      ph(1),shift_phase(1),ph(2),shift_phase(2),phase_alt,phase_alt2);
    disp(string);
  case 3
    disp(string)
    disp('****************************')
    string = sprintf('%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f',...
      ph(1),shift_phase(1),ph(2),shift_phase(2),ph(3),shift_phase(3));
    disp(string);
    
end

%---------------------------------------------------------------------
%
% Plotting
%
function phasingPlot(n_pairs,t,igram,w,abs_array,phase_array,i_max,i_fit);

%ranges for zoom and dots for plotting
i_zoom_f = (i_max-4):(i_max+4);
i_dots = i_fit;%(i_max-1):(i_max+1);

%
%    plot time domain
%
figure(1)
subplot(2,1,1)
plot(t,igram)
axis([t(1) t(end) -1.1 1.1])
set(gca,'XAxisLocation','top')

i_zoom_t = find(t>=-10 & t <= 10);
subplot(2,1,2)
plot(t(i_zoom_t),igram(:,i_zoom_t))
axis([t(i_zoom_t(1)) t(i_zoom_t(end)) -1.1 1.1])

%
%    plot freq domain
%
figure(2),clf
% phase full
h(1)=subplot('position',[.1 .8 .8 .15]);%
plot(w,phase_array)
xlabel('')
ylabel('')
set(h(1),'YTicklabel',[],'Color',[1 1 0.9],...
  'XAxisLocation','top')
% abs full
h(2) = subplot('position',[.1 .5 .8 .3]);
plot(w,abs_array)
xlabel('')
ylabel('')
set(h(2),'XTickLabel',[],'YTicklabel',[])
%phase zoom
h(3)=subplot('position',[.1 .3 .8 .2]);%
plot(w(i_zoom_f),phase_array(:,i_zoom_f))
hold on
hdots = plot(w(i_dots),phase_array(1,i_dots),'o');
set(hdots,'MarkerEdgeColor',[0 0 1],...
  'MarkerFaceColor',[0 0 1])
set(h(3),'YLim',[-pi pi])
if n_pairs>=2
  hdots = plot(w(i_dots),phase_array(2,i_dots),'o');
  set(hdots,'MarkerEdgeColor',[0 .5 0],...
    'MarkerFaceColor',[0 0.5 0])
end
if n_pairs>=3
  hdots = plot(w(i_dots),phase_array(3,i_dots),'o');
  set(hdots,'MarkerEdgeColor',[1 0 0],...
    'MarkerFaceColor',[1 0 0])
end
hold off
xlabel('')
ylabel('')
set(h(3),'XTickLabel',[],'YTicklabel',[],'Color',[1 1 0.9])
%abs zoom
h(4) = subplot('position',[.1 .05 .8 .25]);
plot(w(i_zoom_f),abs_array(:,i_zoom_f))
hold on
hdots = plot(w(i_dots),abs_array(1,i_dots),'o');
set(hdots,'MarkerEdgeColor',[0 0 1],...
  'MarkerFaceColor',[0 0 1])
if n_pairs>=2
  hdots = plot(w(i_dots),abs_array(2,i_dots),'o');
  set(hdots,'MarkerEdgeColor',[0 .5 0],...
    'MarkerFaceColor',[0 0.5 0])
end
if n_pairs>=3
  hdots = plot(w(i_dots),abs_array(3,i_dots),'o');
  set(hdots,'MarkerEdgeColor',[1 0 0],...
    'MarkerFaceColor',[1 0 0])
end
hold off
xlabel('\omega/2\pic (cm-1)')
set(h(4),'YTicklabel',[])

%%compare the abs and real only phased FT
%figure(3)
%plot(w,abs_array(1,:),...
%  w,-real(exp(-1i*phase_array(1,:)).*abs_array(1,:)),...
%  w,-imag(exp(-1i*phase_array(1,:)).*abs_array(1,:)))
