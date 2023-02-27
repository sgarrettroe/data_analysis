function out = fromAbsorptiveToRandNR(w1,w3,absorptive,options)
% go from PP to R and NR spectra and sum

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
ind_mid = floor((ind1(end)+ind1(1))/2);
ind1_ = (1:opt.n_w) - opt.n_w/2 + ind_mid;
w1_ = w1; %this needs to be calculated via dw1 
w3_ = w3; %this needs to be calculated via dw1 

%
[W1,W3] = meshgrid(w1_,w3_);
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

n_t1 = ceil(size(W1,2)/2);
n_t3 = ceil(size(W3,1)/2);

t1 = 0:(size(W1,2)-1);
t3 = 0:(size(W3,1)-1);

w = x; %cropped w1 axis

S1 = sgrsifft2(spec);
S1 = fliplr(circshift(S1,[0 -1]));
S2 = sgrsifft2(spec);

S1(n_t3+1:end,:) = 0;
S1(:,n_t1+1:end) = 0;

S2(:,n_t1+1:end) = 0;
S2(n_t3+1:end,:) = 0;

if opt.flag_plot
    figure,clf,my2dPlot(t1,t3,real(S1),'pumpprobe',false)
    figure,clf,my2dPlot(t1,t3,real(S2),'pumpprobe',false)
end

S1 = sgrsfft2(S1);
S1 = fliplr(circshift(S1,[0 -1]));
S1 = fftshift(S1);

S2 = sgrsfft2(S2);
S2 = fftshift(S2);

if opt.flag_plot
    figure,clf,my2dPlot(W1,W3,real(S1),'pumpprobe',false)
    title('Extracted Re[R_r]')
    
    figure,clf,my2dPlot(W1,W3,real(S2),'pumpprobe',false)
    title('Extracted Re[R_{nr}]')
end

% add for the phase info
S = S1 + S2;

if opt.flag_plot
    figure,clf,my2dPlot(W1,W3,real(S),'pumpprobe',false)
    title('Absorptive spectrum Re[R_{r} + R_{nr}]')
end

out.w1 = w;
out.w3 = w;
out.S = S;
out.R = S1;
out.NR = S2;
