t2_array = [data.t2]; % This takes only the nonempty elements
n_t2_array = length(t2_array);

globalfit.ind1 = find(data(1).w1>globalfit.range1(1) & data(1).w1<globalfit.range1(2));
globalfit.ind3 = find(data(1).w3>globalfit.range3(1) & data(1).w3<globalfit.range3(2));
        
w1 = data(1).w1(globalfit.ind1);
w3 = data(1).w3(globalfit.ind3);
temp = struct;
temp.data = zeros(length(globalfit.ind3),length(globalfit.ind1),n_t2_array); %initialize the variable
temp.sigma = zeros(size(temp.data));
temp.sigma2 = zeros(size(temp.data));

%%
count = 0;
for ii = 1:length(data)
    count = count + 1;
    
    %grab the data in the region of interest
    temp.signal = data(ii).R(globalfit.ind3,globalfit.ind1);
    temp.noise = data(ii).noise(globalfit.ind3,globalfit.ind1);
    
    %normalize to the deepest peak (should be the 0->1 band)
    scale = abs(min(temp.signal(:))); %temp2(:) flattens it to a vector
    
    %save as a slice of the matrix <data>
    temp.data(:,:,count) = (temp.signal)./scale;
    
    %try using the median value of the estimated noise level
    temp.noise = median(temp.noise(:)); %temp3(:) flattens it to a vector
    
    %normalize the noise also
    temp.sigma(:,:,count) = (temp.noise)./scale; % creates a proxy for relative noise based 
    % on the median noise from our more detailed s/n analysis
end
temp.sigma2 = (temp.sigma).^2;

%%
