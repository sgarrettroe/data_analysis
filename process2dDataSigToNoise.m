%% cd to the wherever the data is and load the data
clear data
process2dDataInitialize

cd(localdir{1})

% load designated _summary.mat file
file_name = sprintf('%s-%02i_summary.mat',datestring{1},experiment{1});
load(file_name);
%loading complete

A = exist('data','var');
if ~A
    error('There is no "data" variable.');
end
clear A

%%

phase_steps = 2;
for ii = 1:length(data)
    blargh = sprintf('== Run No. %i ==',ii);
    disp(blargh);
    clear blargh
    datestring2 = data(ii).datestring;
    tempdir = fullfile(datadirbasename,datestring2,'\temp');
    thistempdir = tempdir{1};
    cd(thistempdir)
    nScans = data(ii).PARAMS.nScans;
        if nScans < 4
            error('Your spectrum has too few scans for this script to work.');
        end
    oldR = data(ii).R; %this is the last (cumulative) spectrum

    %set up our array
    dR = zeros([floor(nScans/phase_steps),size(oldR)]); %initializes matrix
    count = 0;
    scanNumber = data(ii).scan_number;
    disp('Remaining scans on this run:');
    for iii = nScans - phase_steps:-phase_steps:1
        count = count + 1;

        m = load(sprintf('%03i-%04i.mat',scanNumber,iii));
        disp(iii)
        %calculate the difference of the accumulated signal from the
        %last phase_steps scans (usually 2). Note that the spectra have
        %to be reweighted for the number of scans that contributed
        dR(count,:,:) = oldR.*(iii + phase_steps) - m.data.R*iii;

        %save the current spectrum for subtracting next time through
        %the loop
        oldR = m.data.R;
    end
        %calculate the standard deviation across all the temp files. Note
        %that the std works on the first non-singleton dimension, which is
        %why we put them in the order dR(count,:,:). The sqrt in the
        %denominator I am not sure about...
    data(ii).noise = squeeze(std(dR))./sqrt(phase_steps);
end
clear oldR m dR count phase_steps nScans scanNumber datestring2 tempdir ...
    thistempdir
disp('finished')
%% Display figures of signal, noise, and signal to noise
for ii = 1:length(data)
    n_contours = 10;
    map=myMapRGB2(n_contours);
      
        figure(1),clf,
        subplot(2,2,1)
        contourf(data(ii).w1(ind1),data(ii).w3(ind3),data(ii).R(ind3,ind1)),title('Signal');
        colormap(map)
        subplot(2,2,2)
        contourf(data(ii).w1(ind1),data(ii).w3(ind3),data(ii).noise(ind3,ind1),12),title('Noise');
        subplot(2,2,3)
        contourf(data(ii).w1(ind1),data(ii).w3(ind3),abs(data(ii).R(ind3,ind1))./data(ii).noise(ind3,ind1),12),title('Signal/Noise') %signal:noise --> abs = absolute value
        pause
end

%%
cd(localdir{1})
beep
flag_save = input('Save output? (1 = "YES"; 0 = "NO"): ');
if flag_save
    save(sprintf('%s-%02i_summary.mat',datestring{1},experiment{1}),'data','-append')
end
flag_save = 0;
beep