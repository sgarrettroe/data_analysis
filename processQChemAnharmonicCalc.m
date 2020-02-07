function out = processQChemAnharmonicCalc(log_name,basis_options,roptions)

flag_output_matrix = 0;

base_name = regexprep(log_name,'.out$',''); %match at the end

%
% read basis options
% 
energy_cutoff = basis_options.energy_cutoff;
nstates = boptions.nstates;
active_modes = boptions.active_modes; % which modes will be treated

r_states = {1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:3,...
    1:3,...
    1:3}; %reduced basis size for each active mode (only define for the active modes)

%
% set up the local modes from the log file and basis options
%

[sod, tod, fod,m,~] = readQChemAnharmonicCalculations(log_name);

count = 0;
n_tod = size(tod,1);
for ii = 1:n_tod
   inds = tod(ii,1:3);
   val = tod(ii,4);
   if inds(1)==inds(2)&inds(1)==inds(3)
       count = count+1;
       tod_diag(count) = val;
   end
end
count = 0;
n_fod = size(fod,1);
for ii = 1:n_fod
   inds = fod(ii,1:4);
   val = fod(ii,5);
   if inds(1)==inds(2)&inds(1)==inds(3)&inds(1)==inds(4)
       count = count+1;
       fod_diag(count) = val;
   end
end

disp('diagonal second, third, and fourth derivatives');
fprintf(1,'%10.4f %10.4f %10.4f\n',[sod;tod_diag;fod_diag]);

%
% assuming inputs of sod, tod, fod (second third and fourth order
% derivatives)
%
% also need mode_reduced_mass (reduced mass list) and dipole directions
% (dmu_dqs = n x 3)


%set up local modes
lmoptions = struct('nstates',[],'m',[],'w',[],'dmu_dq',[]);
n_active_modes = length(active_modes);
count = 0;
for ii = 1:n_active_modes
    count = count+1;
    lmoptions(count).nstates = nstates(active_modes(ii));
    lmoptions(count).m = m(active_modes(ii));
    lmoptions(count).w = sod(active_modes(ii));
    lmoptions(count).dmu_dq(1:3) = dmu_dq(active_modes(ii),1:3);
end

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% rotate to the scf basis
scfmodes = vscf2(lmodes,tod,fod,active_modes);

% reduce the size of the basis for each active mode
rscfmodes = reduceBasisSize(scfmodes,r_states);

% make the product basis
prscfmodes= constructProductModesSparse(rscfmodes);

% scf reduced better matrix elements
pmodes5 = evaluateAnharmonicity2(prscfmodes,tod,fod,active_modes);

rpmodes = reduceProductBasisSize(pmodes5,energy_cutoff);

if flag_output_matrix
    writeMatrixForBarry(base_name);
end

%
% what!?! this is actually crazy sparse at the bottom! Can't I just select
% out these low energy states and mix them!?!
BW = roptions.BW;
zpe = pmodes5.H(1,1);
one_exciton = zpe + sod(6);
two_exciton = zpe + 2*sod(6);

E_harm = sort(diag(pmodes5.H));
n_states = pmodes5.NSTATES;

figure(1),clf
plot(E_harm,'o')
hold on

%g.s.
n_zero_ex = 1;

plot([1 n_states],[zpe zpe])

%one exciton branch
ll = one_exciton - BW/2;
ul = one_exciton + BW/2;
n_one_ex = full(sum(E_harm>=ll & E_harm<=ul));

plot([1 n_states],[ll ll],[1 n_states],[ul ul])

%two exciton branch
ll = two_exciton - 2*BW/2;
ul = two_exciton + 2*BW/2;
n_two_ex = full(sum(E_harm>=ll & E_harm<=ul));

plot([1 n_states],[ll ll],[1 n_states],[ul ul])
hold off

fprintf(1,'n0 \t n1 \t n2\ttot\n');
fprintf(1,'%i\t%i\t%i\t%i\n',n_zero_ex,n_one_ex,n_two_ex,full(sum(E_harm<energy_cutoff)));


%[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

% calculate response
% out(1) = responseFunctions2(pmodes,roptions);
% out(2) = responseFunctions2(pmodes2,roptions);
% out(3) = responseFunctions2(pmodes3,roptions);
% out(4) = responseFunctions2(pmodes4,roptions);
out(1) = responseFunctions2(rpmodes,roptions);
% out(6) = responseFunctions2(pmodes6,roptions);


if flag_plot
% plots
range = rptions.range;
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
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20,'zlimit',0.1);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])
end

    
% difference of vscf energies (E_harm) and CI energies
E_anh = out(1).E;
figure(2),clf
plot(E_harm(1:451)-E_harm(1),'o')
hold on
plot(E_anh,'x')
hold off

end

