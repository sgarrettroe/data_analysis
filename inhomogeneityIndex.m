function [val,val2d] = inhomogeneityIndex(w1,w3,R,NR,peak_pos)
% extract the inhomogeneity index from rephasing and non-rephasing spectra
% see Roberts JCP 2006

%inputs
z_R = abs(R);
z_NR = abs(NR);
x = w1;
y = w3;
%peak_pos = [2.340830964899873   2.335883365565500]*1e3; %w1 w3 position

x0 = peak_pos(1);
y0 = peak_pos(2);
[~,ind_x] = min((x-x0).^2);
[~,ind_y] = min((y-y0).^2);

val2d = (z_R-z_NR)./(z_R+z_NR);
val = val2d(ind_y,ind_x);

