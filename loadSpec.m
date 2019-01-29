%base_name = 'h2o-2';
old_dir = pwd;
switch base_name
  case 'hod_jobfile2'
    cd 'P:/matlab/md/hod/'
    load([base_name,'_spectra.mat']);
    return;
  case 'h2o-2_jobfile3'
    cd 'P:/matlab/md/hod/'
    load([base_name,'_spectra.mat']);
    return;
  case 'hod_jobfile1'
    cd 'c:\Documents and Settings\s.garrett-roe\hod\run-1'
    nt = 32;
    ntint = 4;
    n_steps = 100001;
    time_file_name = cell(1,32);
    n_time_files = 32;
    dt = 0.01;
    for i = 1:n_time_files
      time_file_name{i} = ['times-' num2str(i-1) '.dat'];
    end
  case 'h2o-2'
    cd p:\matlab\md\h2o-2
    dt = 0.01;
    nt  = 32;
    ntint = 6;
    t2_array = [0:6:66, 6:6:108];
    t4_array = [0:6:66 repmat(0,1,18)];
    n_steps = 100001;
  case 'h2o-2_ntint4_'
    cd p:\matlab\md\h2o-2
    dt = 0.01;
    nt  = 32;
    ntint = 4;
    %t2_array = 0:4:12;
    %t4_array = t2_array;
    n_steps = 100001;
    time_file_name ={ ...
      'times-pcp1fa.dat' ...
      'times-pcp1fb.dat' ...
      'times-pcp2fa.dat' ...
      'times-pcp2fb.dat' ...
      'times-pcp1fc.dat' ...
      'times-pcp1fd.dat' ...
      'times-pcp2fc.dat' ...
      'times-pcp2fd.dat' ...
      'times-phononfa.dat' ...
      'times-phononfb.dat' ...
      'times-vibronfc.dat' ...
      'times-vibronfd.dat'}% ...
%      'times-polaronfc.dat' ...
%      'times-polaronfd.dat'};
%16 16 dup	'times-vibronfa.dat' ... 

    %older version
%    time_file_name ={'times-pcp1fa.dat' ...
%      'times-pcp1fb.dat' ...
%      'times-pcp2fa.dat' ...
%      'times-pcp2fb.dat' ...
%      'times-pcp1fb.dat' ...
%      'times-vibronfa.dat'};    
  case 'test_r5'
    cd 'c:\Documents and Settings/s.garrett-roe/md/test_r5/'
    dt = 0.01;
    nt = 32;
    ntint = 6; %how many time points preintegrated
    t2_array = 0;
    t4_array = t2_array;
    n_steps = 1000001;
  case 'test_r5_ntint4'
    cd 'c:\Documents and Settings/s.garrett-roe/md/test_r5/'
    dt = 0.01;
    nt = 32;
    ntint = 4; %how many time points preintegrated
    t2_array = 0:4:20;
    t4_array = t2_array;
    n_steps = 1000001;
    time_file_name = ''
  case 'test_r5_long'
    cd c:/md/test_r5/
    dt = 0.01;
    nt = 64;
    ntint = 6; %how many time points preintegrated
    t2_array = 0;
    t4_array = 0;
    n_steps = 1000001;
  case 'test_r5_fine'
    cd c:/md/test_r5/
    dt = 0.01;
    nt = 32;
    ntint = 3; %how many time points preintegrated
    t2_array = 0;
    t4_array = 0;
    n_steps = 1000001;
  case 'test_r5_coarse'
    cd c:/md/test_r5/
    dt = 0.01;
    nt = 32;
    ntint = 12; %how many time points preintegrated
    t2_array = 0;
    t4_array = 0;
    n_steps = 1000001;
  case 'cum'
    cd p:\matlab\md\test_r5\
    t2_array = 0;
    t4_array = 0;
    n_steps = 0;
    nt = 32;
    dt = .01;
    ntint = 6;
  case 'hod'
    cd p:\matlab\md\temp\
    t2_array = 0;
    t4_array = 0;
    n_steps = 0;
    nt = 32;
    dt = .01;
    ntint = 4;
  otherwise
    error('unknown basename, can''t set parameters nt and ntint')
end

%n_bins = 64; %for joint probability density
%t = (0:99)*dt; %for c2 plotting
w0 = 3450; %wavenumbers -- center frequency

flag_joint_probability = 0;
flag_print = 0;
flag_movie = 0;

%%
%---------------------------------------------------------------------------
%
%    Load spectra
%
%----------------------------------------------------------------------------
%try to load t2_array and t4_array from the times-* files
if exist('time_file_name','var') && ~isempty(time_file_name)
      disp('using time_file_name to load t2 and t4');
  if iscell(time_file_name)
    %if it is a cell array, load each file in order
    t2_t4_pairs = zeros(0,2);
    for i = 1:length(time_file_name)
      t2_t4_pairs = [t2_t4_pairs; load(time_file_name{i})];
    end
  else
    %if not a cell array just load one file
    t2_t4_pairs = load(time_file_name);
  end
  %unpack pairs into t2_ t4_array
  t2_array = t2_t4_pairs(:,1)';
  t4_array = t2_t4_pairs(:,2)';
else
  disp('no time_file_name, looking for t2_array and t4_array...');
  if exist('t2_array','var') ...
      && exist('t4_array','var')...
      && (length(t2_array)==length(t4_array))
    disp('okay, using t2_array and t4_array')
    n_t2_array = length(t2_array);
    t2_t4_pairs = [t2_array(:) t4_array(:)];
  end
end
n_t2_array=length(t2_array);

%these are global variables now
%c=2.99792458e10; %cm/s
%wavenumbersToInvPs=c*1e-12; % ~ 3e-2
%invPsToWavenumbers = 1/wavenumbersToInvPs;% ~ 33

t = (0:nt-1)*(dt*ntint);
freq = fftFreqAxis(t,'time_units','ps','freq_units','cm-1','fftshift','on')
% a=1/(ntint*dt)*invPsToWavenumbers;
% dfreq = a/(nt-1);
% freq = (-a/2:dfreq:a/2+dfreq)-dfreq/2;%+w0; %this is not quite right...
freq = [freq,-freq(1)]; %copy the first (negative) to the positive frequency side

%
% load 1D spectrum
%
spec1d = load([base_name,'_spec1D_t2_0_t4_0.dat']);
spec1d = [spec1d; spec1d(1)];
spec1d = fliplr(spec1d');
%spec1d=flipud(spec1d);

%
% load 2D spectra
%
clear spec2d;
P1=zeros(nt+1,nt+1);
P2=zeros(nt+1,nt+1);
for i=1:n_t2_array;
  %axis
  dummy = load([base_name,'_spec2D_t2_',num2str(t2_array(i)),'_t4_',num2str(t4_array(i)),'.dat']);
  spec2d{i}=zeros(nt,nt);
  count=0;
  for j=1:nt
    for k=1:nt
      count=count+1;
      %spec2d{i}(j,k)=dummy(count,1)+dummy(count,2);
      P1(j,k) = dummy(count,1);
      P2(j+1,k) = dummy(count,2);
    end
  end
  %  P2 = circshift(P2,[1 0]);
  P1(nt+1,:) = P1(1,:);
  P1(:,nt+1) = P1(:,1);
  P2(1,:) = P2(nt+1,:);
  P2(:,nt+1) = P2(:,1);
  spec2d{i} = P1+P2;
  spec2d{i} = flipud(spec2d{i});
%  spec2d{i} = circshift(spec2d{i},[0 -1]);
  spec2d{i} = permute(spec2d{i},[2 1]);
end
clear P1 P2

%%
% load 3D spectra
%
clear spec3d;
R1 = zeros(nt+1,nt+1,nt+1);
R2 = zeros(nt+1,nt+1,nt+1);
R3 = zeros(nt+1,nt+1,nt+1);
R4 = zeros(nt+1,nt+1,nt+1);
for i=1:n_t2_array;
  dummy = load([base_name,'_spec3D_t2_',num2str(t2_array(i)),'_t4_',num2str(t4_array(i)),'.dat']);
  spec3d{i}=zeros(nt+1,nt+1,nt+1);
  count =0;
  for j=1:nt;
    for k=1:nt;
      for l=1:nt;
        count=count+1;
        R1(j,k,l)=dummy(count,1);
        R2(j+1,k,l)=dummy(count,2);
        R3(j,k+1,l)=dummy(count,3);
        R4(j+1,k+1,l)=dummy(count,4);
       %spec3d{i}(j,k,l)=dummy(count,1)+dummy(count,2)+dummy(count,3)+dummy(count,4);
      end
    end
  end
  R1(nt+1,:,:) = R1(1,:,:);
  R1(:,nt+1,:) = R1(:,1,:);
  R1(:,:,nt+1) = R1(:,:,1);
  %
  R2(1,:,:) = R2(nt+1,:,:);
  R2(:,nt+1,:) = R2(:,1,:);
  R2(:,:,nt+1) = R2(:,:,1);
  %
  R3(nt+1,:,:) = R3(1,:,:);
  R3(:,1,:) = R3(:,nt+1,:);
  R3(:,:,nt+1) = R3(:,:,1);
  %
  R4(1,:,:) = R4(nt+1,:,:);
  R4(:,1,:) = R4(:,nt+1,:);
  R4(:,:,nt+1) = R4(:,:,1);
  spec3d{i} = R1 + R2 + R3 + R4;
  %spec3d{i} = circshift(spec3d{i},[0 0 -1]);
  %spec3d{i} = flipdim(spec3d{i},2);
  spec3d{i} = permute(spec3d{i},[2,3,1]);
  %spec3d{i} = circshift(spec3d{i},[-1 -1 -1]);
end

clear R1 R2 R3 R4

% damn it! I still have errors on the extreme points of the spectra (first
%and last), so for now, just truncate everything.
freq=freq(2:end-1);
spec1d =spec1d(2:end-1);
for i = 1:n_t2_array
  spec2d{i}=spec2d{i}(2:end-1,2:end-1);
  spec3d{i}=spec3d{i}(2:end-1,2:end-1,2:end-1);
end

%% go back to the original directory
cd(old_dir);


%% pack them in a structure
s = specStruct;

%% More statistics to extract q_135 from spectra, etc
ss = specStats(freq,spec1d,spec2d,spec3d);

%%
loadSpec_plot1d
loadSpec_plot2d
loadSpec_plot3d

%% slices of 3D
%loadSpec_plot3dSlices
