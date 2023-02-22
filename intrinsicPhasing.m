function phase = intrinsicPhasing(w1,w3,S,options)
% calculate the "intrinsic" phase from complex 2D spectrum
% Philip JM Johnson, Klemens L Koziol, Peter Hamm Optics Express 2017
% https://doi.org/10.1364/OE.25.002928


%opt.peak_pos = [2340.4, 2335.1; 2340.4, 2308.8];
opt.peak_pos = [];
opt.flag_plot = false;
opt.flag_two_level_system = false;
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
x = w1(ind1);
y = w3(ind3);
uh = S(ind3,ind1);

if size(opt.peak_pos,1)==1
    opt.flag_two_level_system = true;
end
if opt.flag_two_level_system
    x0 = opt.peak_pos(1,1);
    y0s = opt.peak_pos(1,2);
else
    x0 = opt.peak_pos(1,1);
    y0s = opt.peak_pos(2,1:2);
end

midpoint = sum(y0s)/2;
[~,i_cut] = min((x-x0).^2);
[~,i_mid] = min((y-midpoint).^2);

if opt.flag_plot
    figure,clf
    plot(y,unwrap(angle(uh(:,i_cut))))
    hold on
    plot([y(1),y(end)],-3*pi/2*[1,1])
    plot([y(1),y(end)],-pi/2*[1,1])
    plot([y(1),y(end)],pi/2*[1,1])
    if opt.flag_two_level_system
        plot(y,-pi*real(uh(:,i_cut))./min(real(uh(:,i_cut))))
    else
        plot([midpoint,midpoint],[-2*pi,2*pi])
        plot(y,pi*real(uh(:,i_cut))./max(real(uh(:,i_cut))))
    end
    hold off
    set(gca,'ylim',[-2*pi,pi])
end

phase = angle(uh(i_mid,i_cut)) + pi/2;
