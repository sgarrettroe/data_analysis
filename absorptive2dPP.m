function [s,PP] = absorptive2dPP(s,varargin)
%calculate the absorptive spectrum from pump-probe 2d data
%
% s = absorptive2d(s,'Property',value,...)
%
% s = absorptive2d(s,'phase',val)
%     uses a phase of val
%
% s = absorptive2d(s,'zeropad',val)
%     The zeropadded length. Should be equal to twice the number of time
%     points for the optimum amount of information in the real spectrum
%     (the default).
%
% s = absorptive2d(s,'range',[lim_l lim_u])
%     Plots over the frequency window of interest given by the lower and
%     upper limits
%
% s = absorptive2d(s,'fft_type','type')
%     Type can be 'fft', the normal fft, or 'sgrsfft' which scales the
%     first data point by 0.5
%
% s = absorptive2d(s,'apodization','type')
%     Can be triangular or gaussian. Others can be implemented by adding
%     the methods to the apodization_list and then changing the window_fxn
%
% s = absorptive2d(s,'pumpprobe',true)
%     The default behavior is to plot 'pump-probe' style
%
% s = absorptive2d(s,'pumpprobe',false)
%     The plots are (x,y) = (omega_1, omega_3)  style
%
% s = absorptive2d(s,'plot',true)
% s = absorptive2d(s,'plot',false)
%     Turn all the plots on or off
%
% s = absorptive2d(s,'plotraw',true)
% s = absorptive2d(s,'plotraw',false)
%     Turn the plots of the raw data on or off
%
% [s,PP] = absorptive2dPP(...) returns the processed time domain data

%default values
n_contours = 20;
zeropad = s.zeropad;
phase = s.phase;
t0_bin = s.t0_bin;
range = [s.freq(1) s.freq(end)];
fft_type = s.fft_type;
fft_type_list = {'fft','sgrsfft'};
apodization = 'none';
apodization_list = {'none','triangular','gaussian'};
flag_pumpprobe = false; %plot style is (w1,w3) as (x,y)
flag_plot=true;
flag_plotraw = false;
flag_fftshift = 'on';
time = s.time;
flag_freq_domain_filter = false;
% freq_domain_filter_center = 2050;
freq_domain_filter_fwhm = 400;

%determine which version of the input arguments are being passed based on
%if the first value is a property string or a phase
if nargin>1
    if isa(varargin{1},'char')
        input_arguments_version = 2;
    else
        input_arguments_version = 1;
    end
    
    switch input_arguments_version
        case 1
            if nargin >= 2
                phase = varargin{1}(1);
                %phase2 = varargin{1}(2);
                s.phase = phase;
                %s.phase2 = phase2;
            end
            if nargin >=3
                if ~isempty(varargin{2})
                    zeropad = varargin{2};
                end
            end
            if nargin >=4
                if ~isempty(varargin{3})
                    range = varargin{3};
                end
            end
            if nargin >=5
                fft_type = varargin{4};
            end
            if nargin >= 6
                if ~isempty(varargin{5})
                    apodization = varargin{5};
                end
            end
            
        case 2
            while length(varargin)>=2
                arg = varargin{1};
                val = varargin{2};
                switch lower(arg)
                    case 'n_contours'
                        n_contours = val;
                        if mod(n_contours,2)
                            warning('SGRLAB:ContourLines','my2dPlot: Odd number of contour lines may produce unexpected results!');
                        end
                    case 'phase'
                        phase = val(1); %take only the first element if it is an array
                    case 'zeropad'
                        zeropad = val;
                    case 'range'
                        range = val;
                    case 'fft_type'
                        fft_type = val;
                    case 'apodization'
                        apodization = val;
                    case {'pumpprobe_style','pumpprobe'}
                        flag_pumpprobe = val;
                    case 'plot'
                        flag_plot = val;
                    case {'freq_domain_filter','filter'}
                        flag_freq_domain_filter = val;
                    case 'filter_center'
                        flag_freq_domain_filter = true;
                        freq_domain_filter_center = val;
                    case 'filter_fwhm'
                        flag_freq_domain_filter = true;
                        freq_domain_filter_fwhm = val;
                    otherwise
                        error(['my2dPlot: unknown option ',arg])
                end
                varargin = varargin(3:end);
            end
    end
end
n_freq = length(s.freq);
if n_freq == 0
    flag_spectrometer = false;
    flag_remove_DC=false;
    flag_plot=false;
else
    flag_spectrometer = true;
    flag_remove_DC=true;
end
n_time = length(s.time);



%error checking of inputs here?
if ~any(strcmpi(fft_type,fft_type_list)), error(['fft type ',fft_type,' not known in absorptive2d.m']);end
if ~any(strcmpi(apodization,apodization_list)), error(['apodization type ',apodization,' not known in absorptive3d.m']);end
s.comment = [s.comment,' fft_type ',fft_type,' apodization ',apodization];

if flag_spectrometer
    %begin calculation
    
    switch apodization
        case 'none'
            window_fxn = ones(1, n_time);
        case 'triangular'
            window_fxn = linspace(1,0,n_time);
        case 'gaussian'
            window_fxn = exp(-(linspace(0,3,n_time)).^2);
    end
    window_fxn = repmat(window_fxn,n_freq,1);
    
    w1 = fftFreqAxis(time,'time_units',s.time_units,...
        'freq_units',s.freq_units,...
        'fftshift',flag_fftshift,...
        'zeropad',zeropad);
    R = zeros(zeropad,n_freq);
    PP = s.PP;
    
    if flag_remove_DC
        ncol = size(PP,2);
        PP = PP - repmat(mean(PP,2),1,ncol);
    end
    
    %all data before 0 gets set to 0 because it is from a mixture of
    %different population times
    preface = PP(:,1:t0_bin-1);
    PP(:,1:t0_bin-1) = 0;
    
    %make it the mean rather than 0 to reduce zeropadding artifacts. The
    %alternative would be to subtract the mean from all elements of the PP
    %data.
    %PP(:,1:t0_bin-1) = repmat(mean(PP(:,t0_bin:end),2),1,t0_bin-1);
    
    %rotate bins to put t0 as the first bin
    PP = circshift(PP,[0 -t0_bin+1]);
    if flag_freq_domain_filter
        
        w1_ = fftFreqAxis(time,'time_units',s.time_units,...
        'freq_units',s.freq_units,...
        'fftshift','on','zeropad',2*size(PP,2));
    
        nrows = length(s.w3);
        
        %make a fourier filter (a top hat)
        filter_center = mean(s.freq);
%         filter_ind = (w1_ >= filter_center - freq_domain_filter_fwhm & w1_ <= filter_center + freq_domain_filter_fwhm);
        filter_ind = (w1_ >= filter_center - freq_domain_filter_fwhm & w1_ <= filter_center + freq_domain_filter_fwhm | ...
            w1_ <= -filter_center+freq_domain_filter_fwhm & w1_ >= -filter_center - freq_domain_filter_fwhm);
        filter_fxn = zeros(size(w1_));
        filter_fxn(filter_ind) = 1;
        
        WIN = repmat(filter_fxn,nrows,1);
        
        PP_ = bsxfun(@minus, PP, mean(PP,2)); % subtract the mean to avoid a spike in the FT at 0 frequency
%         R = sgrsfft(PP_,[],2); % have to FT along the correct dimension
        R = fftshift(sgrsfft(PP_,2*size(PP_,2),2),2);
        R = R.*WIN;

%         PP = real(sgrsifft(R,[],2));
        PP = real(sgrsifft(fftshift(R,2),[],2));
        PP = PP(:,1:end/2);
%         PP(:,end-t0_bin+2:end) = 0;
    end
    
    PP_ = PP.*window_fxn;
    switch fft_type
        case 'fft'
            R = fft(PP_,zeropad,2);
        case 'sgrsfft'
            R = sgrsfft(PP_,zeropad,2);
    end
    
    %correct for the phase calculated by the phasing routine
    R = real(R.*exp(-1i*phase*pi/180));
    R = fftshift(R,2);
    
    %   if flag_remove_DC
    %     R(:,1)=0;
    %   end
    
    
    %redo frequency axis in case we zeropadded
    %s = freq2d(s,zeropad);
    %  s = freq2d(s,'zeropad',zeropad,...
    %    'spectrometer',flag_spectrometer,...
    %    'fftshift',flag_fftshift);
    
    %store the results
    s.w1 = w1;
    s.w3 = s.freq;
    s.R = R;
    PP = circshift(PP,[0 t0_bin-1]);
    if ~flag_freq_domain_filter
        PP(:,1:t0_bin-1) = preface;
    end
    s.PP = PP;
    %filter in frequency domain if requested
    %     if flag_freq_domain_filter
    %
    %         nrows = length(s.w3);
    %
    %         %make a window function (a gaussian) to surpress frequencies outside of the
    %         %region of interest. The factor of 2.355 converts fwhm to standard
    %         %deviation
    %         window_fxn = exp(-(w1-freq_domain_filter_center).^2./(2*(freq_domain_filter_fwhm/2.355)^2));
    %
    %         %copy the window_fxn so it is the same size as the data
    %         WIN = repmat(window_fxn,nrows,1);
    %         %calc new PP data by
    %         %1) multiply by window
    %         %2) shift so zero frequency is the first element of the array
    %         %3) inverse fft
    %         %4) take only the real part
    %         PP = real(ifft(ifftshift(s.R.*WIN,2),[],2));
    %
    %         %fix the time axis to match the fft length
    %         dt = time(2) - time(1);
    %         time = (0:(zeropad-1))*dt;
    %     end
    
    %
    %    Plot results
    %
    map = myMapRGB2(n_contours+1);
    ind = find(s.w1>range(1) & s.w1<range(2));
    
    if flag_plot
        %input raw data
        if flag_plotraw
            figure(100),clf
            if flag_pumpprobe
                x = s.freq;
                y = time;
                %      z = real(s.PP)';
                z = real(PP)';
                x_label = '\omega_{probe} / 2\pic';
                y_label = 't_{pump} / fs';
            else
                x = time;
                y = s.freq;
                %      z = real(s.PP);
                z = real(PP);
                x_label = 't_1 / fs';
                y_label = '\omega_3 / 2\pic';
            end
            %    size(x)
            %    size(y)
            %    size(z)
            contourf(x,y,z,n_contours)
            %axis square
            myCaxis2(z,n_contours);
            colormap(map)
            xlabel(x_label)
            ylabel(y_label)
            drawnow
        end
        
        %fft data the "real" spectrum
%         figure(101),clf
        figure(101),clf
        if flag_pumpprobe
            x = s.w3;
            y = s.w1(ind);
            z = s.R(:,ind)';
        else
            x = s.w1(ind);
            y = s.w3;
            z = s.R(:,ind);
        end
        a=my2dPlot(x,y,z,'n_contours',n_contours,'pumpprobe',flag_pumpprobe);
    end %if flag_plot
    
else
    %if time domain experiment/simulation
    error('SGRLAB:notimplemented','absorptive2dPP all time domain analysis not implemented');
end
