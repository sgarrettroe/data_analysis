%on a mac, this should be placed in  ~/Documents/MATLAB
%in windows, this should be .../Documents and Settings/MATLAB

%use this for the path to where the data analysis functions are placed
%addpath ~/Developer/data_analysis

global c c_SI c_cm c_cmfs wavenumbersToInvFs wavenumbersToInvPs ...
  wavenumbersToInvSec invSecToWavenumbers invFsToWavenumbers ...
  invPsToWavenumbers fringeToFs fringeToPs...
  q h hbar epsilon_0 mu_0 debyeToCm N_A ...
  k_B amuToKg
c_SI = 2.9979e8;
c = 2.9979e10;
c_cm = c;
c_cmfs = 2.9979e-5;
wavenumbersToInvFs = c*1e-15;
wavenumbersToInvPs=c*1e-12;
wavenumbersToInvSec = c;
invFsToWavenumbers=1/wavenumbersToInvFs;
invPsToWavenumbers=1/wavenumbersToInvPs;
invSecToWavenumbers = 1/wavenumbersToInvSec;
fringeToFs = 632.8e-7/c_cmfs; %FIXME!!!
fringeToPs = fringeToFs/1000;
q = 1.6e-19; %charge coulombs
h = 6.626e-34; %plancks
hbar = 1.054571596e-34;
epsilon_0 = 8.8541878e-12; %F/m
mu_0 = 1.2566470e-6; %H/m
debyeToCm = 3.33564e-30; %Cm (coulomb meters)
N_A = 6.022e23; 
k_B = 1.38e-12; %J K-1
amuToKg= 1.660538e-27; %kg

