%% Generate toy data

n = 20; %subjects
r = 5; %ROIs
t = 50; %time points
s = 2; %sessions

rng(5);

roi_data = rand(n,r,t,s);
behavior = randi([50,90],1,n); %I used the score range from the CBCL


save('roi_data.mat','behavior')
save('behavior.mat','behavior') 

%% Temporal Pairwise

pairwise_temporal_ISC = zeros(n,n,r); %create matrix (n, n, r)
for i = 1:n
    for j = 1:n
        for roi = 1:r
            dati = roi_data(i,roi,:,:); %flatten into vector i 
            dati = reshape(dati,[1,size(dati,3)*size(dati,4)]);
            datj = roi_data(j,roi,:,:); %flatten into vector j
            datj = reshape(datj,[1,size(datj,3)*size(datj,4)]); 
            pairwise_temporal_ISC(i,j,roi) = corr(dati',datj'); %transpose, correlate i and j, store in output matrix
        end
    end
end

%% Temporal LOO

%zeros are built into function
loo_temporal_ISC = get_loo_isc(roi_data); %uses fuction built at the bottom of the script


%% Dynamic 

w = 5; %5 time windows
step_size = t/w; %50 time points/5 time windows equals window size 10
loo_dynamic_ISC = zeros(n,r,w); %create empty (n,r,w) matrix
for T = 1:w
    t1 = T*step_size - step_size + 1; %select sliding window
    t2 = T*step_size;
    roi_data_T = roi_data(:,:,t1:t2,:); %calculate ISC within each window
    loo_dynamic_ISC(:,:,T) = get_loo_isc(roi_data_T); %store in output matrix
end

%% create image

dyn_isc_im = squeeze(loo_dynamic_ISC(1,1,:)); %remove dimesions of 1
dyn_isc_im = plot(dyn_isc_im); %store plot
saveas(dyn_isc_im,'dynamic_isc_plot.pdf'); %save to Week6 folder

%% Spatial

loo_spatial_ISC = zeros(n,1); %create empty (n,1) matrix
for subs = 1:n
    roi_sub = roi_data(subs,:,:,:); %separate one sub, find mean of the rest
    roi_loo = mean(roi_data(1:size(roi_data,1) ~= subs,:,:,:),1);
    dat_sub = mean(roi_sub,3); %average over time
    dat_sub = mean(dat_sub,4); %average over sessions
    dat_sub = squeeze(dat_sub); %flatten into vector
    dat_loo = mean(roi_loo,3); %same thing for one subject
    dat_loo = mean(dat_loo,4);
    dat_loo = squeeze(dat_loo);
    loo_spatial_ISC(subs) = corr(dat_sub',dat_loo'); %correlate between odd one out and the rest, store in output matrix
end




%% Intra-Subject Correlation

intrasubject_temporal_ISC = zeros(n,r); %create empty (n,r) matrix
for subs = 1:n
    for roi = 1:r
    session1 = squeeze(roi_data(subs,roi,:,1)); %flatten session 1 and 2 into seperate flat vectors
    session2 = squeeze(roi_data(subs,roi,:,2));
    intrasubject_temporal_ISC(subs,roi) = corr(session1,session2); %correlate between sessions, store in output matrix
    end
end


%% ISFC

loo_ISFC = zeros(n,r,r); %create empty (n,r,r) matrix
for subs = 1:n
    roi_loo = mean(roi_data(1:n ~= subs,:,:,:),1); 
    for roi1 = 1:r
        for roi2 = 1:r
            datsubs_roi1 = roi_data(subs,roi1,:,:); %average over all subjects except 1
            datsubs_roi1 = reshape(datsubs_roi1,[t*s,1]);
            dat_loo_roi2 = roi_loo(1,roi2,:,:);
            dat_loo_roi2 = reshape(dat_loo_roi2,[t*s,1]); % flatten into vector
            loo_ISFC(subs,roi1,roi2) = corr(datsubs_roi1,dat_loo_roi2); % correlate and store
        end
    end
end


%% ISRSP to Behavior

behavioral_similarity = zeros(n,n); %create empty (n,n) matrix
%difference between behavioral scores
for i = 1:n
    for j = 1:n
        behavioral_similarity(i,j) = abs(behavior(i) - behavior(j));
    end
end

behavioral_similarity = -behavioral_similarity; % take inverse

results = zeros(r,2); %create resutls matrix with room for r and p for each ROI
for roi = 1:r
    neural_similarity = pairwise_temporal_ISC(:,:,roi);
    n_sim = extractRDM(neural_similarity); %use Lab 4 function `extractRDM`
    b_sim = extractRDM(behavioral_similarity);
    [rho,p] = corr(n_sim,b_sim,'type','Spearman');
    results(roi,1) = rho;
    results(roi,2) = p;
end

%% Save Matrices

save('output.mat','pairwise_temporal_ISC','loo_temporal_ISC','loo_dynamic_ISC','loo_spatial_ISC','intrasubject_temporal_ISC','loo_ISFC','behavioral_similarity','results'); 

%% Calculate LOO ISC

function lootisc = get_loo_isc(data)
n = size(data,1); 
r = size(data,2);
lootisc = zeros(n,r); %create empty matrix
for i = 1:n
    roi_loo = mean(data(1:n~= i,:,:,:),1);
    for roi = 1:r
        dati = data(i,roi,:,:); %average over all subs except one
        dati = reshape(dati,[size(dati,3)*size(dati,4),1]); %flatten
        dat_loo = roi_loo(1,roi,:,:); %repeat for odd one out
        dat_loo = reshape(dat_loo,[size(dat_loo,3)*size(dat_loo,4),1]);
        lootisc(i,roi) = corr(dati,dat_loo); %store correlation
    end
end
end





