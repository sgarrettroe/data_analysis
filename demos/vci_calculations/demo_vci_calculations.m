%% Notes for the group meeting
% the object here is to walk through the ideas behind the coupled
% oscillator calculations

%% Demo of local mode basis
% mode 1
lmoptions(1).nstates = 3; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 10; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

%% the wavefunction in the local mode basis
% here is the first basis vector (ie the ground state)
lmodes(1).u(1)

% here is a way to access the states by v=1 notation (not the ground
% state!)
lmodes(1).v(1)


%% Quiz 1: creation operators

% creation operator answer to quiz 1
lmodes(1).c

%% annihilation
lmodes(1).a

%% position (q) & momentum (p) ops
lmodes(1).q
lmodes(1).p

%%
% show hamiltonian / energies of mode 1. Note that the energies have the zero point
% energies included
edit('constructLocalModes.m')
lmodes(1).h
lmodes(1).h_

%% Demo: Kronecker product is a way to represent a tensor product
v1 = [1; 2; 3; 4]
v2 = [5; 6; 7]

v1v2 = kron(v1,v2)

%% 
a = 1;
b = 0;
c = 0;
e = 0;
f = 0;
g = 1;

v1 = [a; b; c]
v2 = [e; f; g]

v1v2 = kron(v1,v2)


%%
a = sqrt(2)/2;
v1v2 = [0 a 0 a 0 0 0 0 0]'

displayCoeffMatrix(v1v2,lmodes) 

%%

kron([0; 1; 0],[1; 0; 0])

kron([1; 0; 0],[0; 1; 0])

v1v2 = kron([0; 1; 0],[1; 0; 0]) + kron([1; 0; 0],[0; 1; 0])

%% alternate way to write that
v1v2 = 1/sqrt(2).*(kron(lmodes(1).v(1),lmodes(2).v(0))...
    + kron(lmodes(1).v(0),lmodes(2).v(1)))

%% Demo of product modes
% mode 1
lmoptions(1).nstates = 3; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 10; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 3; %number of states in basis for mode 1
lmoptions(2).m = 1; %mass
lmoptions(2).w = 8; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% make the product basis
pmodes = constructProductModes(lmodes);

%upack the results
f=fieldnames(pmodes)
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii}))
end
clear f     

%% Demo of creation and annihilation ops
PSI = pmodes.IDENTITY(:,1)

% action of creation ops
C1*PSI

C2*PSI

C*PSI

%% Couplings between local modes

delta3 = -0; %this adds coupling but very little energy shift. weird
delta4 = -6;

anh_cubic_1 = -0;
anh_cubic_2 = -0;
anh_quartic_1 = -0.87;
anh_quartic_2 = -0.92;

% mode 1
lmoptions(1).nstates = 20; %number of states in basis for mode 1
lmoptions(1).m = 16; %reduced mass?
lmoptions(1).w = 1889.184+10.5+12; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(1).dmu_dq = [sqrt(0.908) 0 0]; %transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 20; %number of states in basis for mode 2
lmoptions(2).m = 16; %reduced mass?
lmoptions(2).w = 1947.621+11.1+12; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [0 sqrt(0.671) 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);


%
% make the product basis
pmodes = constructProductModes(lmodes);
%upack the results
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end
clear f     

%
  H_ = H + anh_cubic_1*(A1+C1)^3 + anh_cubic_2*(A2+C2)^3 ...
      + anh_quartic_1*(A1+C1)^4 + anh_quartic_2*(A2+C2)^4 ...
  + delta3*((A1+C1)^2*(A2+C2) + (A2+C2)^2*(A1+C1)) ...
  + delta4*((A1+C1)^2*(A2+C2)^2);

pmodes.H_ = H_;

[V,E] = analyzeEnergyLevels(lmodes,pmodes,'ind',1:6);


%% try to reproduce Rhcomplex with responseFunctions2
clear lmoptions lmodes pmodes roptions

% one can play with these values and see the spectrum change

% coupling between modes
delta3 = -0; 
delta4 = -6;
%delta3 = -20; 
%delta4 = -2;

% local mode anharmonicities
anh_cubic_1 = -0.8;
anh_cubic_2 = -0.8;
anh_quartic_1 = -1;
anh_quartic_2 = -1;

% mode 1
lmoptions(1).nstates = 10;  % number of states in basis for mode 1
lmoptions(1).m = 16;  % reduced mass
lmoptions(1).w = 1889.184+10.5+12;  % omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(1).dmu_dq = [1 0 0];  % transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 10;  % number of states in basis for mode 2
lmoptions(2).m = 16;  % reduced mass?
lmoptions(2).w = 1947.621+11.1+12;  % omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [-1 0 0];  % transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);


%
% make the product basis
pmodes = constructProductModes(lmodes);

% unpack the results (for demo purposes; not a good idea for production)
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end
clear f     

%
  H_ = H + anh_cubic_1*(A1+C1)^3 + anh_cubic_2*(A2+C2)^3 ...
      + anh_quartic_1*(A1+C1)^4 + anh_quartic_2*(A2+C2)^4 ...
  + delta3*((A1+C1)^2*(A2+C2) + (A2+C2)^2*(A1+C1)) ...
  + delta4*((A1+C1)^2*(A2+C2)^2);

pmodes.H_ = H_;

%[V,E] = analyzeEnergyLevels(lmodes,pmodes,'ind',1:6);
%
%
target = sort([0 1889.184 1947.621 3767.517 3883.911 3813.698]');

%
% setup response
roptions.order = 3;
roptions.energy_cutoff = [];
roptions.flag_energy_cutoff = false;
roptions.w0 = 1915;
roptions.n_t = 128;
roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.dt = 0.25;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 1915;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',2);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);
fprintf(1,'\n\n')
fprintf(1,'%f \t%f \t%f\n',[target E(1:6)-E(1) (-target + E(1:6)-E(1))]')

%
% calculate response
out = responseFunctions2(pmodes,roptions);

%
range = [1850 1980];
ind1 = (out.w1 >= range(1) & out.w1 <= range(2));
ind3 = (out.w3 >= range(1) & out.w3 <= range(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind1,ind3);
J = out.J(ind1);
figure(101),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])

%% try to reproduce figures from the book about fermi resonance A) W/O anharmonicity REDUX
% compare to Edler JCP 2003 Fig 3a
% essential features are there -- inversion of the esa peaks, single
% "cross-peak". Nevertheless, my peaks are more complex in shape. Are there
% pathways (like coherence transfer peaks) that I am including that they
% did not that cause a phase twisted appearance? 
%
% ok interesting. Turning off the transitions from the bend (setting mu to
% 0) gives exactly the shapes in the paper...
clear lmoptions lmodes pmodes roptions

% mode 
lmoptions(1).nstates = 10; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 1760/2; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 10; %number of states in basis for mode 2
lmoptions(2).m = 1; %mass
lmoptions(2).w = 1760; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% experimental parameters
BW = 300; %bandwidth
w_laser = 1760; %central frequency

% setup response
roptions.order = 3;
roptions.w0 = 1760;
%roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = w_laser;
roptions.BW = BW;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;


% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% make the product basis
pmodes = constructProductModes(lmodes);
%upack the results
f=fieldnames(pmodes)
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'));
end
clear f     

Delta = -15;

%  H_ = H + anh_cubic_1*Q1^3 + anh_cubic_2*Q2^3 ...
%      + anh_quartic_1*Q1^4 + anh_quartic_2*Q2^4

%try dunham style
%pmodes.H_ = H - 1*(C1*A1)^2 - 1*(C2*A2)^2 -Delta*(C1*C1*A2+ A1*A1*C2);
pmodes.H_ = H + Delta*(C1*C1*A2+ A1*A1*C2);

%

[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);

%
% calculate response
tic
out = responseFunctions2(pmodes,roptions);
toc

%
range = [1690 1830];
ind1 = (out.w1 >= range(1) & out.w1 <= range(2));
ind3 = (out.w3 >= range(1) & out.w3 <= range(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind1,ind3);
J = out.J(ind1);
figure(101),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])

%% try to reproduce figures about fermi resonance B) W/ anharmonicity
% compare to Edler JCP 2003 Fig 3b. 
clear lmoptions lmodes pmodes roptions

chi = -.8;
Delta = -15;

% mode 1
lmoptions(1).nstates = 10; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 1760/2; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [0 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 10; %number of states in basis for mode 2
lmoptions(2).m = 1; %mass
lmoptions(2).w = 1770; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% make the product basis
pmodes = constructProductModes(lmodes);
%upack the results
f=fieldnames(pmodes)
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'));
end
clear f     


%  H_ = H + anh_cubic_1*Q1^3 + anh_cubic_2*Q2^3 ...
%      + anh_quartic_1*Q1^4 + anh_quartic_2*Q2^4

%try dunham style
H_ = H - 0*(C1*A1)^2 + chi*(C2+A2)^4 -Delta*(C1*C1*A2+ A1*A1*C2);
%H_ = H + Delta*(C1*C1*A2+ A1*A1*C2);
%H_ = H +0*(C1*A1)^2 - chi*(C2*A2)^2 -Delta*(Q1*Q1*Q2+Q2*Q2*Q1);
pmodes.H_ = H_;

% setup response
roptions.order = 3;
roptions.w0 = 1760;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 1760;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);

[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% calculate response
tic
out = responseFunctions2(pmodes,roptions);
toc

%
range = [1690 1830];
ind1 = (out.w1 >= range(1) & out.w1 <= range(2));
ind3 = (out.w3 >= range(1) & out.w3 <= range(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind1,ind3);
J = out.J(ind1);
figure(101),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])

%% try single inhomogeneously broadened mode
% seems to work out ok
clear lmoptions lmodes pmodes roptions

chi = -2.5;

% mode 1
lmoptions(1).nstates = 20; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 2150; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [0 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% make the product basis
pmodes = constructProductModes(lmodes);
%upack the results
f=fieldnames(pmodes)
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'));
end
clear f     


%  H_ = H + anh_cubic_1*Q1^3 + anh_cubic_2*Q2^3 ...
%      + anh_quartic_1*Q1^4 + anh_quartic_2*Q2^4

%try dunham style
H_ = H + chi*(C1+A1)^4;
pmodes.H_ = H_;

% setup response
roptions.order = 3;
roptions.w0 = 2150;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 50; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2150;
roptions.BW = 300;
roptions.c2form = '2exp1fast';
roptions.c2params = struct('T2',1,'Delta1_cm',5,'Delta2_cm',10,'tau1',2,'tau2',30);

[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

% calculate response
tic
out = responseFunctions2(pmodes,roptions)
toc

%
range = [2025 2200];
ind1 = (out.w1 >= range(1) & out.w1 <= range(2));
ind3 = (out.w3 >= range(1) & out.w3 <= range(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind1,ind3);
J = out.J(ind1);
figure(101),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])
%my2dPlot(w(ind),w(ind),R(ind,ind)','n_contours',20)


%% try to reproduce Clinton's BH4- data
clear lmoptions lmodes pmodes roptions

% chi = -2.5;
% Delta1 = -22;
% Delta2 = -22;

% chi = -1.3;%-2.5;
% Delta1 = -22;
% Delta2 = -28; %-22;

chi = -2.5%-3.043;
Delta1 = -22;
Delta2 = -25;
% mode 1
lmoptions(1).nstates = 10; %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 1090; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [0 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2
lmoptions(2).nstates = 10; %number of states in basis for mode 2
lmoptions(2).m = 1; %mass
lmoptions(2).w = 1190; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 3
lmoptions(3).nstates = 10; %number of states in basis for mode 2
lmoptions(3).m = 1; %mass
lmoptions(3).w = 2250; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(3).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

% make the product basis
pmodes = constructProductModes(lmodes);
%upack the results
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'));
end
clear f     


%  H_ = H + anh_cubic_1*Q1^3 + anh_cubic_2*Q2^3 ...
%      + anh_quartic_1*Q1^4 + anh_quartic_2*Q2^4

%try dunham style
H_ = H - 0*(C1*A1)^3 + chi*(C3+A3)^4 ...
    -Delta1*(C1*C1*A3+ A1*A1*C3)...
    -Delta2*(C1*C2*A3+ A1*A2*C3);
%H_ = H + Delta*(C1*C1*A2+ A1*A1*C2);
%H_ = H +0*(C1*A1)^2 - chi*(C2*A2)^2 -Delta*(Q1*Q1*Q2+Q2*Q2*Q1);
pmodes.H_ = H_;

% setup response
roptions.order = 3;
roptions.w0 = 2225;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0.5; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2225;
roptions.BW = 300;
roptions.c2form = '3exp1fast';
roptions.c2params = struct('T2',4,'Delta1_cm',12,'Delta2_cm',18,'tau1',0.2,'tau2',0.4,'Delta3_cm',8,'tau3',30);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);

% calculate response
tic
out = responseFunctions2(pmodes,roptions);
toc

%%
range = [2130 2320];
ind1 = (out.w1 >= range(1) & out.w1 <= range(2));
ind3 = (out.w3 >= range(1) & out.w3 <= range(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind1,ind3);
J = out.J(ind1);
figure(101),clf
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])

%%
eigenstate_of_interest = 7;
vec = V(:,eigenstate_of_interest);
vec = conj(vec).*vec;
[~,ind] = sort(vec,'descend');
vec = vec(ind);
ind = ind(vec>=0.1);
vec = vec(vec>=0.1);

vs = indexToVs(ind,lmodes);

%set up some output
nmodes = length(lmodes);

formatstring = '(';
for ii = 1:nmodes
    formatstring = strcat(formatstring,'%3d,');
end
formatstring(end) = ')'

fprintf(1,'EIGENSTATE %3d',eigenstate_of_interest)
fprintf(1,'index\tcontribution\tquantum nums\n')
fprintf(1,strcat('%d\t%f\t',formatstring,'\n'),[ind,vec,vs(:,1),vs(:,2),vs(:,3)]')

%% 

%% CO2 input from qchem
logname = 'CO2.out';
[sod, tod, fod,m,dmu_dq] = readQChemAnharmonicCalculations(logname);
nmodes = length(sod);
tod_diag = zeros(1,nmodes);
fod_diag = zeros(1,nmodes);
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

%%
% assuming inputs of sod, tod, fod (second third and fourth order
% derivatives)
%
% also need mode_reduced_mass (reduced mass list) and dipole directions
% (dmu_dqs = n x 3)
clear lmoptions lmodes pmodes roptions

%calculation parameters
nstates = 6; %10 takes ~60 sec per diagonalization

active_modes = [1:4];

%set up local modes
n_active_modes = length(active_modes);

lmoptions = struct('nstates',[],'m',[],'w',[],'dmu_dq',[]);

count = 0;
for ii = 1:n_active_modes
    count = count+1;
    lmoptions(count).nstates = nstates;
    lmoptions(count).m = m(active_modes(ii));%mode_reduced_mass(active_modes(ii)); %FIX FIX FIX
    lmoptions(count).w = sod(active_modes(ii));
    lmoptions(count).dmu_dq(1:3) = dmu_dq(active_modes(ii),1:3);%dmu_dqs(active_modes(ii),:); %FIX FIX FIX
end

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

prod([lmodes.nstates])

% make the product basis
pmodes = constructProductModesSparse(lmodes);
%pmodes = constructProductModes(lmodes);

pmodes = evaluateAnharmonicity(pmodes,tod,fod,active_modes);


%


% experimental parameters
BW = 300; %bandwidth
w_laser = 2380; %central frequency

% setup response parameters
roptions.order = 3;
roptions.w0 = 2450;
roptions.n_t = 1024;
%roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = w_laser;
roptions.BW = BW;
roptions.c2form = '1exp1fast';
roptions.c2params = struct('T2',3,'Delta1_cm',2,'tau1',40);
roptions.T = 300;
roptions.verbose = 1;

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);

%
% calculate response
tic
out = responseFunctions2(pmodes,roptions);
toc

%%
%range = [(w_laser-BW/2) (w_laser+BW/2)];
range1 = [2360 2420];
range3 = [2330 2420];
ind1 = (out.w1 >= range1(1) & out.w1 <= range1(2));
ind3 = (out.w3 >= range3(1) & out.w3 <= range3(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind3,ind1);
J = out.J(ind1);
figure(101),clf
%[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20,'zlimit',0.1);
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])


%% CO2 basis set size dependence
%
% code here down has been copied from previous presentations but not yet
% tested
%
%%
% assuming inputs of sod, tod, fod (second third and fourth order
% derivatives)
%
% also need mode_reduced_mass (reduced mass list) and dipole directions
% (dmu_dqs = n x 3)
clear lmoptions lmodes pmodes roptions

%calculation parameters
nstates = 4; %10 takes ~60 sec per diagonalization

active_modes = [1:4];

%set up local modes
n_active_modes = length(active_modes);

lmoptions = struct('nstates',[],'m',[],'w',[],'dmu_dq',[]);

count = 0;
for ii = 1:n_active_modes
    count = count+1;
    lmoptions(count).nstates = nstates;
    lmoptions(count).m = m(active_modes(ii));%mode_reduced_mass(active_modes(ii)); %FIX FIX FIX
    lmoptions(count).w = sod(active_modes(ii));
    lmoptions(count).dmu_dq(1:3) = dmu_dq(active_modes(ii),1:3);%dmu_dqs(active_modes(ii),:); %FIX FIX FIX
end

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

prod([lmodes.nstates])

% make the product basis
pmodes = constructProductModesSparse(lmodes);
%pmodes = constructProductModes(lmodes);

pmodes = evaluateAnharmonicity(pmodes,tod,fod,active_modes);


%


% experimental parameters
BW = 300; %bandwidth
w_laser = 2380; %central frequency

% setup response parameters
roptions.order = 3;
roptions.w0 = 2450;
roptions.n_t = 1024;
%roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = w_laser;
roptions.BW = BW;
roptions.c2form = '1exp1fast';
roptions.c2params = struct('T2',3,'Delta1_cm',2,'tau1',40);
roptions.T = [];
roptions.verbose = 1;

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);

%
% calculate response
tic
out = responseFunctions2(pmodes,roptions);
toc

%
range = [(w_laser-BW/2) (w_laser+BW/2)];range1=range;range3=range;
%range1 = [2360 2420];
%range3 = [2330 2420];
ind1 = (out.w1 >= range1(1) & out.w1 <= range1(2));
ind3 = (out.w3 >= range3(1) & out.w3 <= range3(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind3,ind1);
J = out.J(ind1);
figure(101),clf
%[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20,'zlimit',0.1);
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])

beep

%% Examples of matrix elements
q_op = lmodes(1).q;
home
disp('X^1')
q_op
lmodes(1).q_pow_1

disp('X^2')
q_op*q_op
lmodes(1).q_pow_2

disp('X^3')
q_op*q_op*q_op
lmodes(1).q_pow_3

disp('X^4')
q_op*q_op*q_op*q_op
lmodes(1).q_pow_4


%% BETTER MATRIX ELEMENTS

clear lmoptions lmodes pmodes* roptions out

D111 = 100; %cubic
D112 = -200;
D221 = -200;
D222 = 100;
D1111 = -100;%quartic
D2222 = -300;

%build tod and fod
tod = [1 1 1 D111;...
    1 1 2 D112;...
    2 2 2 D222
    2 2 1 D221];

fod = [1 1 1 1 D1111;
    2 2 2 2 D2222];

n_states = [12 12];
r_states = {[1:5], 1:3};
active_modes = 1:2;

% mode 1
lmoptions(1).nstates =n_states(1); %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 2350/2; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [0 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2 
lmoptions(2).nstates = n_states(2); %number of states in basis for mode 2
lmoptions(2).m = 1; %mass
lmoptions(2).w = 2350; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates



% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);
rlmodes = reduceBasisSize(lmodes,r_states);

scfmodes = vscf2(lmodes,tod,fod);
rscfmodes = reduceBasisSize(scfmodes,r_states);

%figure(1),plot(diag(D)-diag(scfmodes(i_mode).h))


%
% make the product basis
pmodes = constructProductModesSparse(lmodes);
prmodes= constructProductModesSparse(rlmodes);
pscfmodes= constructProductModesSparse(scfmodes);
prscfmodes= constructProductModesSparse(rscfmodes);

% full
pmodes     = evaluateAnharmonicity(pmodes,tod,fod,active_modes);
% reduced with better matrix elements
pmodes2    = evaluateAnharmonicity2(prmodes,tod,fod,active_modes);
% reduced with worse matrix elements
pmodes3    = evaluateAnharmonicity(prmodes,tod,fod,active_modes);
% % scf full
% pmodes4  = evaluateAnharmonicity(pscfmodes,tod,fod,active_modes);
% % scf reduced better
% pmodes5 = evaluateAnharmonicity2(prscfmodes,tod,fod,active_modes);
% % scf reduced worse
% pmodes6 = evaluateAnharmonicity(prscfmodes,tod,fod,active_modes);

% setup response
roptions.order = 3;
roptions.w0 = 2250;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2250;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);

%[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

% calculate response
out(1) = responseFunctions2(pmodes,roptions);
out(2) = responseFunctions2(pmodes2,roptions);
out(3) = responseFunctions2(pmodes3,roptions);
% out(4) = responseFunctions2(pmodes4,roptions);
% out(5) = responseFunctions2(pmodes5,roptions);
% out(6) = responseFunctions2(pmodes6,roptions);

%
range = [2100 2500];
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

%% MOTIVATE VSCF
% CO2 input from qchem
logname = 'CO2.out';
[sod, tod, fod,m,dmu_dq] = readQChemAnharmonicCalculations(logname);
nmodes = length(sod);
tod_diag = zeros(1,nmodes);
fod_diag = zeros(1,nmodes);
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
clear lmoptions lmodes pmodes roptions

%calculation parameters
nstates = 6; %10 takes ~60 sec per diagonalization

active_modes = [1:4];

%set up local modes
n_active_modes = length(active_modes);

lmoptions = struct('nstates',[],'m',[],'w',[],'dmu_dq',[]);

count = 0;
for ii = 1:n_active_modes
    count = count+1;
    lmoptions(count).nstates = nstates;
    lmoptions(count).m = m(active_modes(ii));%mode_reduced_mass(active_modes(ii)); %FIX FIX FIX
    lmoptions(count).w = sod(active_modes(ii));
    lmoptions(count).dmu_dq(1:3) = dmu_dq(active_modes(ii),1:3);%dmu_dqs(active_modes(ii),:); %FIX FIX FIX
end

% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);

prod([lmodes.nstates])

% make the product basis
pmodes = constructProductModesSparse(lmodes);
%pmodes = constructProductModes(lmodes);

pmodes = evaluateAnharmonicity(pmodes,tod,fod,active_modes);

%%
% the harmonic hamiltonian
lmodes(4).h
e_harm = diag(lmodes(4).h);

% the anharmonic hamiltonian
h_ = lmodes(4).h + 1/6*tod_diag(4)*lmodes(4).q_pow_3 + 1/24*fod_diag(4)*lmodes(4).q_pow_4

[v,e] = eig(h_,'vector')

n_states = lmodes(4).nstates;
x_plot = 1:n_states;
figure(1),clf
subplot(2,1,1)
plot(x_plot,e_harm,'bo-',x_plot,e,'rx-')
subplot(2,1,2)
plot(x_plot,e-e_harm,'r-')

% note that the 'diagonal' anharmonicity gets the difference of the 0-1 and
% 1-2 spacing with the wrong sign!!!

%% CONTINUE VSCF IDEAS 

clear lmoptions lmodes pmodes* roptions out

D111 = 100; %cubic
D112 = -200;
D221 = -200;
D222 = 100;
D1111 = -100;%quartic
D2222 = -300;

%build tod and fod
tod = [1 1 1 D111;...
    1 1 2 D112;...
    2 2 2 D222
    2 2 1 D221];

fod = [1 1 1 1 D1111;
    2 2 2 2 D2222];

n_states = [12 12];
r_states = {[1:5], 1:3};
active_modes = 1:2;

% mode 1
lmoptions(1).nstates =n_states(1); %number of states in basis for mode 1
lmoptions(1).m = 1; %mass
lmoptions(1).w = 2350/2; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
% lmoptions(1).dmu_dq = [0 0 0]; %transition dipole in cartesian coordinates
lmoptions(1).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates

% mode 2 
lmoptions(2).nstates = n_states(2); %number of states in basis for mode 2
lmoptions(2).m = 1; %mass
lmoptions(2).w = 2350; %omega/(2 pi c) in wavenumbers (cm-1) (gets converted to angular freq in units of rad/ps)
lmoptions(2).dmu_dq = [1 0 0]; %transition dipole in cartesian coordinates



% make our local mode basis sets
lmodes = constructLocalModes(lmoptions);
rlmodes = reduceBasisSize(lmodes,r_states);

scfmodes = vscf2(lmodes,tod,fod);
rscfmodes = reduceBasisSize(scfmodes,r_states);


%
% make the product basis
pmodes = constructProductModesSparse(lmodes);
prmodes= constructProductModesSparse(rlmodes);
pscfmodes= constructProductModesSparse(scfmodes);
prscfmodes= constructProductModesSparse(rscfmodes);

% full
pmodes     = evaluateAnharmonicity(pmodes,tod,fod,active_modes);
% reduced with better matrix elements
pmodes2    = evaluateAnharmonicity2(prmodes,tod,fod,active_modes);
% reduced with worse matrix elements
pmodes3    = evaluateAnharmonicity(prmodes,tod,fod,active_modes);
% scf full
pmodes4  = evaluateAnharmonicity(pscfmodes,tod,fod,active_modes);
% scf reduced better
pmodes5 = evaluateAnharmonicity2(prscfmodes,tod,fod,active_modes);
% scf reduced worse
pmodes6 = evaluateAnharmonicity(prscfmodes,tod,fod,active_modes);

% setup response
roptions.order = 3;
roptions.w0 = 2250;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2250;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);

%[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

% calculate response
out(1) = responseFunctions2(pmodes,roptions);
out(2) = responseFunctions2(pmodes2,roptions);
out(3) = responseFunctions2(pmodes3,roptions);
out(4) = responseFunctions2(pmodes4,roptions);
out(5) = responseFunctions2(pmodes5,roptions);
out(6) = responseFunctions2(pmodes6,roptions);

%
range = [2100 2500];
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

%% Looks like we have a recipe! do the bh4 hamiltonian

%
cd ~/'Box Sync'/Projects/borohydride/quantum_chemistry

% from qchem
logname = 'bh4-.out';
[sod, tod, fod,m,dmu_dq] = readQChemAnharmonicCalculations(logname);

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

%%
% assuming inputs of sod, tod, fod (second third and fourth order
% derivatives)
%
% also need mode_reduced_mass (reduced mass list) and dipole directions
% (dmu_dqs = n x 3)
clear lmoptions lmodes pmodes* roptions 

flag_output_matrix = 0;

%calculation parameters
nstates = [12 12 12 12 12 12 12 12 12];
active_modes = [1:8];
r_states = {1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:3,...
    1:3,...
    1:3}; %only define for the active modes


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

count = 0;
for ii = 1:n_active_modes
    count = count+1;
    lmoptions(count).nstates = nstates(ii);
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
%pmodes = constructProductModesSparse(lmodes);
%prmodes= constructProductModesSparse(rlmodes);
%pscfmodes= constructProductModesSparse(scfmodes);
prscfmodes= constructProductModesSparse(rscfmodes);

%
% full
%pmodes     = evaluateAnharmonicity(pmodes,tod,fod,active_modes);
% reduced with better matrix elements
%pmodes2    = evaluateAnharmonicity2(prmodes,tod,fod,active_modes);
% reduced with worse matrix elements
%pmodes3    = evaluateAnharmonicity(prmodes,tod,fod,active_modes);
% scf full
%pmodes4  = evaluateAnharmonicity(pscfmodes,tod,fod,active_modes);
% scf reduced better matrix elements
pmodes5 = evaluateAnharmonicity2(prscfmodes,tod,fod,active_modes);
% scf reduced worse
%pmodes6 = evaluateAnharmonicity(prscfmodes,tod,fod,active_modes);

%output pmodes.h_
if flag_output_matrix
    [ii,jj,S]=find(pmodes5.H_);
    string = sprintf('bh4-_vscf_hamiltonian_matrix.csv');
    FID = fopen(string,'w+');
    for kk = 1:length(ii)
        fprintf(FID,'%d,%d,%d\n',ii(kk),jj(kk),S(kk));
    end
    fclose(FID);
end

fprintf(1,'NSTATES = %i\n',pmodes5.NSTATES)

if pmodes5.NSTATES>5000
    fprintf(1,'Number of states %i probably too large...\n',pmodes5.NSTATES);
    return
end
%
% setup response
roptions.order = 3;
roptions.w0 = 2250;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2250;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);

%[V,E]=analyzeEnergyLevels(lmodes,pmodes);

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

% calculate response
% out(1) = responseFunctions2(pmodes,roptions);
% out(2) = responseFunctions2(pmodes2,roptions);
% out(3) = responseFunctions2(pmodes3,roptions);
% out(4) = responseFunctions2(pmodes4,roptions);
out(1) = responseFunctions2(pmodes5,roptions);
% out(6) = responseFunctions2(pmodes6,roptions);

%
range = [2100 2500];
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

%% FROM IMPORT BARRY'S RESULT
%% 
cd ~/'Box Sync'/Projects/borohydride/quantum_chemistry

%% from qchem
logname = 'bh4-.out';
[sod, tod, fod,m,dmu_dq] = readQChemAnharmonicCalculations(logname);

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
clear lmoptions lmodes pmodes* roptions out

flag_output_matrix = 0;

%calculation parameters
nstates = [12 12 12 12 12 12 12 12 12];
active_modes = [1:8];
r_states = {1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:5,...
    1:3,...
    1:3,...
    1:3}; %only define for the active modes


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
%pmodes = constructProductModesSparse(lmodes);
%prmodes= constructProductModesSparse(rlmodes);
%pscfmodes= constructProductModesSparse(scfmodes);
prscfmodes= constructProductModesSparse(rscfmodes);

%
% full
%pmodes     = evaluateAnharmonicity(pmodes,tod,fod,active_modes);
% reduced with better matrix elements
%pmodes2    = evaluateAnharmonicity2(prmodes,tod,fod,active_modes);
% reduced with worse matrix elements
%pmodes3    = evaluateAnharmonicity(prmodes,tod,fod,active_modes);
% scf full
%pmodes4  = evaluateAnharmonicity(pscfmodes,tod,fod,active_modes);
% scf reduced better matrix elements
pmodes5 = evaluateAnharmonicity2(prscfmodes,tod,fod,active_modes);
% scf reduced worse
%pmodes6 = evaluateAnharmonicity(prscfmodes,tod,fod,active_modes);

%% 
% output pmodes5.H_ matrix to file. Send to Barry. 
%% read in the eigenvals and vecs from Barry's PETSc / MPI code (full 85k x 85k matrix)


eigenvals = csvread('eigenvalues_sgr.txt');%_sgr trims the first line of text off so just numerical vals
n_eigenvals = length(eigenvals);

eigenvectors = csvread('eigenvectors.csv');

%
eigenvectors = eigenvectors';
n_eigenvectors = size(eigenvectors,2);
l_eigenvector = size(eigenvectors,1);
if n_eigenvals ~= n_eigenvectors
    error('number of eigenvals %i does not match number of eigenvectors %i',n_eigenvals,n_eigenvectors);
end


%% plot the results
E_harm = diag(pmodes5.H);
energy_cutoff = 120000;

BW = 500;
zpe = eigenvals(1);
one_exciton = zpe + sod(6)
two_exciton = zpe + 2*sod(6)

% = sort(diag(eigenvals));
n_states = n_eigenvals;

figure(1),clf
plot(eigenvals,'o')
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

%%
pmodesBarry = pmodes5;
pmodesBarry.V = eigenvectors;
pmodesBarry.E = eigenvals;

%
clear out roptions
% setup response
roptions.order = 3;
roptions.w0 = 2250;
% roptions.n_t = 1024;
roptions.n_t = 128;
%roptions.n_zp = 2*roptions.n_t; % the zeropadded length
roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[1; 0; 0],[1; 0; 0],[1;0;0],[1;0;0] };
roptions.w_laser = 2250;
roptions.BW = 300;
roptions.c2form = '1fast';
roptions.c2params = struct('T2',0.75);


% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

out = responseFunctions3(pmodesBarry,roptions);

%
%range = [2150 2400];
range = [2100 2500];
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

