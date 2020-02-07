function [out,pmodesBarry] = processBarrysResult( fname_eigenvalues,fname_eigenvectors,pmodes)
%PROCESSBARRYSRESULT Read and plot spectrum from full matrix
%diagonalization

[V,E] = readBarrysResult(fname_eigenvalues,fname_eigenvectors);

pmodesBarry = pmodes;
pmodesBarry.V = V; %potential transpose problem here !!! @@@
pmodesBarry.E = E; 

out = responseFunctions3(pmodesBarry,roptions);

%
%range = [2150 2400];
range = roptions.range;
w1 = out(1).w1;
w3 = out(1).w3;
ind1 = (w1 >= range(1) & w1 <= range(2));
ind3 = (w3 >= range(1) & w3 <= range(2));
w1 = w1(ind1);
w3 = w3(ind3);
for ii = 1:length(out)
R = out(ii).R(ind1,ind3);
J = out(ii).J(ind1);
figure(100+ii),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])
end

end

function [ eigenvectors,eigenvalues ] = readBarrysResult( fname_eigenvalues,fname_eigenvectors )
%READBARRYSRESULT Read the files from barry
%   Detailed explanation goes here

row_offset = 1; %skip the first row
eigenvalues = csvread(fname_eigenvalues,row_offset);%_sgr trims the first line of text off so just numerical vals
n_eigenvals = length(eigenvalues);

eigenvectors = csvread(fname_eigenvectors);

%
eigenvectors = eigenvectors';
n_eigenvectors = size(eigenvectors,2);
l_eigenvector = size(eigenvectors,1);
if n_eigenvals ~= n_eigenvectors
    error('number of eigenvals %i does not match number of eigenvectors %i',n_eigenvals,n_eigenvectors);
end


end

