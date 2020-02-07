function [aniso,iso] = calc2DIRAnisotropy(para,perp,varargin)

%defaults
range_para = [min(para(1).w3) max(para(1).w3)];
range_perp = [min(perp(1).w3) max(perp(1).w3)];

while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'range_para'
            range_para = val;
        case 'range_perp'
            range_perp = val;
    end
    varargin = varargin(3:end);
end

aniso = para;
iso = para;

n_spectra = length(para);

for ii = 1:n_spectra
    
    indpar  = find(para(ii).w3<range_para(2) & para(ii).w3>range_para(1));
    indperp = find(perp(ii).w3<range_perp(2) & perp(ii).w3>range_perp(1));
    
    iso(ii).R = para(ii).R(indpar,:) + 2*perp(ii).R(indperp,:);
    aniso(ii).R = (para(ii).R(indpar,:) - perp(ii).R(indperp,:))./iso(ii).R;
    
    iso(ii).w3 = para(ii).w3(indpar);
    aniso(ii).w3 = para(ii).w3(indpar);
    
end

