function varargout = inhomogeneityIndex(varargin)
% function [val,val2d] = inhomogeneityIndex(w1,w3,R,NR,peak_pos)
% extract the inhomogeneity index from rephasing and non-rephasing spectra
% see Roberts JCP 2006
n_args = 4;
USE_PEAK = 1;
USE_ROI = 2;

if length(varargin)<3
    error(['Not enough input arguments. Call as',...
        'result = inhomogeneityIndex(w1,w3,absorptive,options);',...
        ' or ',...
        'result = inhomogeneityIndex(d.w1,d.w3,d.R,options);'])
end
peak_pos = varargin{end};
varargin = varargin(1:end-1);
n_spectra = length(varargin)/n_args;

varargin = reshape(varargin,n_spectra,n_args)';

switch length(peak_pos)
    case 2
        flag_use_case = USE_PEAK;
    case 4
        flag_use_case = USE_ROI;
    otherwise
        error('expected peak_pos input argument to be either length 2 or 4');
end
if nargout == 1
    flag_struct_out = true;
else
    flag_struct_out = false;
end

if flag_struct_out
    out(n_spectra) = struct('val',[],'w1',[],'w3',[],'val2d',[],'ffcf',[],'ffcf2d',[]);
else
    % if array out
    val_out = zeros(1,n_spectra);
    val2d_out = repmat(zeros(size(varargin{3})),[1,1,n_spectra]);
    ffcf_out  = zeros(1,n_spectra);
end

for i_loop = 1:n_spectra
    
    w1 = varargin{1};
    w3 = varargin{2};
    R = varargin{3};
    NR = varargin{4};

    %inputs
    z_R = abs(R);
    z_NR = abs(NR);
    x = w1;
    y = w3;
    %peak_pos = [2.340830964899873   2.335883365565500]*1e3; %w1 w3 position
    

    val2d = (z_R-z_NR)./(z_R+z_NR);

    if flag_use_case==USE_PEAK
        x0 = peak_pos(1);
        y0 = peak_pos(2);
        [~,ind_x] = min((x-x0).^2);
        [~,ind_y] = min((y-y0).^2);
    end
    if flag_use_case==USE_ROI
        x0 = peak_pos(1);
        y0 = peak_pos(2);
        x1 = peak_pos(3);
        y1 = peak_pos(4);
        ind_x = x >= x0 & x<=x1;
        ind_y = y >= y0 & y<=y1;
        [ind_x,ind_y] = meshgrid(ind_x,ind_y);
    end
    
    val = sum(val2d(ind_y,ind_x),'all');
    
    if flag_struct_out
        out(i_loop).val = val;
        out(i_loop).ffcf = sin(pi/2*val);
        out(i_loop).val2d = val2d;
        out(i_loop).ffcf2d = sin(pi/2*val2d);
        out(i_loop).w1 = w1;
        out(i_loop).w3 = w3;
    else
        val_out(i_loop) = val;
        val2d_out(:,:,i_loop) = val2d;
    end
    
    varargin = varargin(:,2:end);
end

if flag_struct_out
    varargout{1} = out;
else
    varargout{1} = val_out;
    varargout{2} = val2d_out;
end

