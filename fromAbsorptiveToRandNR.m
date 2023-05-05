function out = fromAbsorptiveToRandNR(varargin)
%function out = fromAbsorptiveToRandNR(w1,w3,absorptive,options)
% go from PP to R and NR spectra and sum

if length(varargin)<3
    error(['wrong number of input parameters call as',
        'result = fromAbsorptiveToRandNR(w1,w3,absorptive,options);',...
        ' or ',
        'result = fromAbsorptiveToRandNR(d.w1,d.w3,d.R,options);'])
end
options = varargin{end};
varargin = varargin(1:end-1);
n_spectra = length(varargin)/3;

varargin = reshape(varargin,n_spectra,3)';

out(n_spectra).w1 = [];
out(n_spectra).w3 = [];
out(n_spectra).S = [];
out(n_spectra).R = [];
out(n_spectra).NR = [];
out(n_spectra).t2 = [];

for i_loop = 1:n_spectra
    w1 = varargin{1};
    w3 = varargin{2};
    absorptive = varargin{3};
    
    % defaults
    opt.n_w = 256; %number of freq points
    opt.phase = 0;
    opt.flag_plot = false;
    opt.range1 = [w3(1) w3(end)];
    opt.range3 = [w3(1) w3(end)];
    
    % assign properties from the input options struct
    props = fields(opt);
    for ii = 1:length(props)
        if isfield(options,props{ii})
            opt.(props{ii}) = options.(props{ii});
        end
    end
    
    ind1 = find(w1>=opt.range1(1) & w1<=opt.range1(2));
    ind3 = find(w3>=opt.range3(1) & w3<=opt.range3(2));
    
    % put the spectrum in a bigger freq space in w3
    % is this needed???
    %ind_mid = floor((ind1(end)+ind1(1))/2);
    %dw = w1(2)-w1(1);
    %ind1_ = (1:opt.n_w) - opt.n_w/2 + ind_mid;
    %ind1_ = ind1_(ind1_>0&ind1_<=length(ind1));
    %w1_ = ind1_*dw + w1(ind_mid);
    
    %
    %[W1,W3] = meshgrid(w1_,w1_);
    %x = w1(ind1_);
    [W1,W3] = meshgrid(w1,w1);
    x = w1(ind1);
    y = w3(ind3);
    z = real(absorptive(ind3,ind1)*exp(-1i*opt.phase));
    big_spec = interp2(x,y,z,W1,W3,'makima',0);
    
    
    if opt.flag_plot
        figure
        my2dPlot(W1,W3,big_spec,'pumpprobe',false)
        title('big spectrum')
    end
    
    %
    spec = ifftshift(big_spec);
    
    %n_t = opt.n_w/2;
    n_t = floor(length(x)/2);
    t = 0:(2*n_t-1);
    w = x;
    
    S1 = sgrsifft2(spec);
    S1 = fliplr(circshift(S1,[0 -1]));
    S2 = sgrsifft2(spec);
    
    S1(n_t+1:end,:) = 0;
    S1(:,n_t+1:end) = 0;
    
    S2(:,n_t+1:end) = 0;
    S2(n_t+1:end,:) = 0;
    
    if opt.flag_plot
        figure,clf,my2dPlot(t,t,real(S1),'pumpprobe',false)
        figure,clf,my2dPlot(t,t,real(S2),'pumpprobe',false)
    end
    
    S1 = sgrsfft2(S1);
    S1 = fliplr(circshift(S1,[0 -1]));
    S1 = fftshift(S1);
    
    S2 = sgrsfft2(S2);
    S2 = fftshift(S2);
    
    if opt.flag_plot
        figure,clf,my2dPlot(w,w,real(S1),'pumpprobe',false)
        title('Extracted Re[R_r]')
        
        figure,clf,my2dPlot(w,w,real(S2),'pumpprobe',false)
        title('Extracted Re[R_{nr}]')
    end
    
    % add for the phase info
    S = S1 + S2;
    
    if opt.flag_plot
        figure,clf,my2dPlot(w,w,real(S),'pumpprobe',false)
        title('Absorptive spectrum Re[R_{r} + R_{nr}]')
    end
    
    out(i_loop).w1 = w;
    out(i_loop).w3 = w;
    out(i_loop).S = S;
    out(i_loop).R = S1;
    out(i_loop).NR = S2;
    
    varargin = varargin(:,2:end);
end
