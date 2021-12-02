%% Feature Extraction algorithm to extract psd based features from finger flex ECoG data.

%% Parameters
sub_name = 'bp';
freq_band_1 = "70_135";  % "4-8","8-12","12-40","40-70","70-135","135-200"
freq_band_2 = "135_200";  % "4-8","8-12","12-40","40-70","70-135","135-200"

freq_band_list = ["4_8","8_12","12_40","40_70","70_135","135_200"];
low_cutoff_psd = [2,3,4,11,19,35]; 
upper_cutoff_psd = [3,4,11,18,34,51];  

Index_1= find(contains(freq_band_list,freq_band_1));
low_cut_1 = low_cutoff_psd(Index_1)
upper_cut_1 = upper_cutoff_psd(Index_1)

Index_2= find(contains(freq_band_list,freq_band_2));
low_cut_2 = low_cutoff_psd(Index_2)
upper_cut_2 = upper_cutoff_psd(Index_2)


%% Load data

% Load dataset .mat file
ECoG_Handpose = load(strcat(".\datasets\", sub_name, "\", sub_name, "_fingerflex.mat")); 
data = ECoG_Handpose.data;      %ECoG recordings
Fs = 1000;                      %Sampling Frequency
stim = ECoG_Handpose.cue;       %Paradiagm Signal 
num_channels = size(data,2)     %Obtain number of recorded ECoG channels

%% Common average re-referencing
ca_filtered = com_avg(data,num_channels);  % Common average re-referencing

% Obtain frequency domain representation and visualize
ca_fft = fft(ca_filtered(:,2));
[n,p] = size(ca_filtered(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;

figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Common Average Filtered ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% 6th order cascaded notch filtering
% Sixth order cascaded Butterworth notch filters were applied at frequencies 60 Hz, 120 Hz and 180 Hz.

notch_filtered=ca_filtered;
[n,m]=size(notch_filtered);
for i = 1:6
    Fc=[60*i-0.2,60*i+0.2];     %Notch filter at 60 hz and its harmonics
    Fc
    [b,a] = butter(3,Fc/(Fs/2),'stop');

    %freqz(b,a);
    for j =1:m
        notch_filtered(:,j)=filtfilt(b,a,notch_filtered(:,j)) ;

    end
end
 
% Obtain frequency domain representation and visualize
notch_fft=fft(notch_filtered(:,num_channels));
[n,p]= size(notch_filtered(:,:));
P2 = abs(notch_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Notch filtered ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Spectral Whitening 
% Data is whitened by applying the spectral whitening filter proposed by Grunewald et al., 2019 
%to equalize the spectral contributions.

for i =1:m        
    notch_filtered(:,i)= spectral_whitening(notch_filtered,10,i) ;
end

% Obtain frequency domain representation and visualize
notch_fft=fft(notch_filtered(:,num_channels));
[n,p]= size(notch_filtered(:,:));
P2 = abs(notch_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;

figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Notch filtered ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Apply band pass filter for all channels  (0.5-200)
% The data is bandpass filtered (Chebyshev type II, 3rd order) 
% between 0:5 Hz and 200 Hz.

[n,m]=size(notch_filtered);
for i =1:m
     notch_filtered(:,i)=high_pass(notch_filtered(:,i),0.5,Fs) ;
     notch_filtered(:,i)=low_pass(notch_filtered(:,i),210,Fs);            
end


ca_fft=fft(notch_filtered(:,num_channels));
[n,p]= size(notch_filtered(:,2));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Bandpass Filtered ECoG (After Whitening)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% plotting paradigm and dataglove data (obtain derivative and peak detection on glove signal)

[n,p]= size(data);
t=1:n;
glove_signal = ECoG_Handpose.flex;   % Load glove signals. 
stim_norm = (stim - min(stim))./(max(stim) - min(stim));
glove_norm = (glove_signal - min(glove_signal,[],1))./(max(glove_signal,[],1) - min(glove_signal,[],1));
derv_glove = diff(glove_norm);
% derv_glove=abs(derv_glove);
peaks_glove=[];
peaks=[];

%Peak detection on the derivation of glove signals to detect onset.
for i=1:5
    
    [qrspeaks_X,locs_X] = findpeaks(derv_glove(:,i),t(2:end),'MinPeakHeight',0.01,...               %bp= 0.04  ca=0.03
        'MinPeakDistance',1);
    temp=zeros(10000,1);
    temp(1:length(locs_X),1)=locs_X-1;
    peaks_glove=[ peaks_glove temp];
    
    temp2=zeros(10000,1);
    temp2(1:length(qrspeaks_X),1)=qrspeaks_X;
    peaks=[ peaks temp2];
end

%Visualization
figure;
plot(t,stim_norm);  
hold on;
plot(t,glove_norm)
legend("paradigm",'Thumb','Index','Middle','Ring','Little')
xlabel('Time'), ylabel('data');
title('Finger glove signal with paradigm signal');

figure;
plot(t,stim_norm);  
hold on;
plot(t(2:end),derv_glove);
legend("paradigm",'Thumb','Index','Middle','Ring','Little')
xlabel('Time'), ylabel('data');
title('Derivative of finger glove signal with paradigm signal');

figure;
plot(t,stim_norm);  
hold on;
plot(peaks_glove,peaks,'+');
legend("paradigm",'Thumb','Index','Middle','Ring','Little')
xlabel('Time'), ylabel('data');
title('Detected peaks with paradigm signal');


%% Segment ECoG segments for Thumb

ecog_data=notch_filtered();
relax_times=extracttask(stim,1);   %%% thumb==1
sum(relax_times)/2000
[n,m]=size(ecog_data);
thumb_data=[];                     %change here
data_num=0;
i=1;
while i<=n-1900
    if relax_times(i)~=0
        k=i+299;   %
        first=0;
        thumb =0;
        index=0;
        middle=0;
        ring=0;
        little=0;
        while k<=i+1699                         % checking whther the subject has done the task correctly or not using glove signal
            if ismember(k,peaks_glove(:,1))     %for thumb            
                thumb=1;
                if first==0
                    start=k;
                    first=1;
                end
                break
            end
            if ismember(k,peaks_glove(:,2))      %for index finger
                index=1;
            end
            if ismember(k,peaks_glove(:,3))     %for middle finger
                middle=1;
            end
            if ismember(k,peaks_glove(:,4))     %for ring finger
                ring=1;
            end
            if ismember(k,peaks_glove(:,5))     %for little finger
                little=1;
            end
            k=k+1;
        end
        if thumb==1 && index==0 && middle==0 && ring ==0&& little ==0   %for Thumb
            data_num=data_num+1;
            disp(["Time Segment No (Thumb):" num2str(data_num)]);
            pause(2);
            a=start-2000;
            thumb_data(:,:,data_num)=ecog_data(a:a+3999,:);           
        end
        i=i+2099 ;
    end
    i=i+1;
    
end

%% Temporal feature extraction for Thumb
temporal_thumb=[];

n=size(thumb_data,3);
for i=1:n
    disp(["Feature Extraction for Time Segment No :" num2str(i)])
    pause(2);
    channel=[];
    for j=1:num_channels                 %for each channel
        channel_task=[];
        channel_relax=[];
        
        start=1;
        for k=1:76                    %take periodogram for relax time
            window=thumb_data(start:start+249,j,i); %250ms window    
            [Pxx_relax,F_relax] = periodogram(window,hamming(length(window)),length(window),Fs);
            start=start+50;   %stride 50 ms
            Pxx_relax=10*log10(Pxx_relax);
            channel(j,k)=mean(Pxx_relax(low_cut_1:upper_cut_1));        %(2:18) 2:3 3:4 4:11 11:18 19:34 35:51
            channel(j+num_channels,k)=mean(Pxx_relax(low_cut_2:upper_cut_2));   %(19:76)
            
        end
        
    end
    
    temporal_thumb(:,:,i)=channel;
    
end


%% data of Thumb

channel_mean=mean(temporal_thumb(:,1:16,:),2);
temporal_thumb_normalized=temporal_thumb-channel_mean;
temporal_thumb_normalized_filt= imfilter(temporal_thumb_normalized,fspecial('gaussian',[3,3],1));
temporal_thumb_normalized_all=mean(temporal_thumb_normalized_filt,3);
temporal_thumb_m=mean(temporal_thumb,3);
temporal_thumb_m= imfilter(temporal_thumb_m,fspecial('gaussian',[1,1],1));

figure;
f=1:num_channels;
t=-2:0.05:2;
imagesc(t,f,temporal_thumb_m(1:num_channels,:));  
grid on
colormap('jet')
title("Feature Map for Thumb in 70-135 Hz Before Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')
figure;
f=1:num_channels;
t=-2:0.05:2;
imagesc(t,f,temporal_thumb_m(num_channels+1:end,:,1));   
grid on
colormap('jet')
title("Feature Map for Thumb in 135-200 Hz Before Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')

figure;
imagesc(t,f,temporal_thumb_normalized_all(1:num_channels,:));
grid on
colormap('jet')
title("after nomalizing Thumb 4-8 Hz data name : wm")
figure;
imagesc(t,f,temporal_thumb_normalized_all(num_channels+1:num_channels*2,:));
grid on
colormap('jet')
title("after nomalizing Thumb 8-12 Hz data name : wm")

%% for Index

ecog_data=notch_filtered();
relax_times=extracttask(stim,2);   %%% Index==2
sum(relax_times)/2000
[n,m]=size(ecog_data);
index_data=[];                     %change here
data_num=0;
i=1;
while i<=n-1900
    if relax_times(i)~=0
        k=i+299;        %start checking from +100ms to 1800ms
        first=0;
        thumb =0;
        index=0;
        middle=0;
        ring=0;
        little=0;
        while k<=i+1699                          % checking whther the subject has done the task correctly or not using glove signal
            if ismember(k,peaks_glove(:,2))      %for index finger
                index=1;
                if first==0
                    start=k;
                    first=1;
                end
                break
            end
            if ismember(k,peaks_glove(:,1))     %for thumb            
                thumb=1;
            end
            
            if ismember(k,peaks_glove(:,3))     %for middle finger
                middle=1;
            end
            if ismember(k,peaks_glove(:,4))     %for ring finger
                ring=1;
            end
            if ismember(k,peaks_glove(:,5))     %for little finger
                little=1;
            end
            k=k+1;
        end
        if  index==1% && thumb==0  && middle==0 && ring ==0 && little ==0  %for index
            data_num = data_num+1;
            disp(["Time Segment No (Index):" num2str(data_num)]);
            pause(2);
            a=start-2000;
            index_data(:,:,data_num)=ecog_data(a:a+3999,:);           
        end
        i=i+2099 ;
    end
    i=i+1;
    
end

%% Temporal feature extraction for index
temporal_index=[];

n=size(index_data,3);
for i=1:n
    disp(["Feature Extraction for Time Segment No (Index):" num2str(i)]);
    pause(2);
    channel=[];
    for j=1:num_channels                 %for each channel
        channel_task=[];
        channel_relax=[];
        
        start=1;
        for k=1:76                    %take periodogram for relax time
            window=index_data(start:start+249,j,i); %250ms window    %%%%%%only change here
            [Pxx_relax,F_relax] = periodogram(window,hamming(length(window)),length(window),Fs);
            start=start+50;   %stride 50 ms
            Pxx_relax=10*log10(Pxx_relax);
            channel(j,k)=mean(Pxx_relax(low_cut_1:upper_cut_1));        %(2:18)  2:3 3:4 4:11 11:18 19:34 35:51
            channel(j+num_channels,k)=mean(Pxx_relax(low_cut_2:upper_cut_2));   %(19:76)
            
        end
        
    end
    
    temporal_index(:,:,i)=channel;
    
end


%% data of Index

channel_mean=mean(temporal_index(:,1:16,:),2);
temporal_index_normalized=temporal_index-channel_mean;
temporal_index_normalized_filt= imfilter(temporal_index_normalized,fspecial('gaussian',[3,3],1));
temporal_index_normalized_all=mean(temporal_index_normalized_filt,3);



figure;
f=1:num_channels;
t=0.125:0.05:3.875;
imagesc(t,f,temporal_index(1:num_channels,:,1));   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("before nomalizing Index Finger")
figure;
imagesc(t,f,temporal_index_normalized_all(1:num_channels,:));
grid on
colormap('jet')
title("after nomalizing Index Finger 4-8 Hz data name : wm")
figure;
imagesc(t,f,temporal_index_normalized_all(num_channels+1:num_channels*2,:));
grid on
colormap('jet')
title("after nomalizing Index Finger 8-12 Hz data name : wm")

%% for Middle

ecog_data=notch_filtered();
relax_times=extracttask(stim,3);   %%% Middle==3
sum(relax_times)/2000
[n,m]=size(ecog_data);
middle_data=[];                     %change here
data_num=0;
i=1;
while i<=n-1900
    if relax_times(i)~=0
        k=i+299;
        first=0;
        thumb =0;
        index=0;
        middle=0;
        ring=0;
        little=0;
        while k<=i+1699%2499%1899        %1899  %2499->jp                 % checking whther the subject has done the task correctly or not using glove signal
            if ismember(k,peaks_glove(:,3))     %for middle finger
                middle=1;
                if first==0
                    start=k;
                    first=1;
                end
                break
            end
            if ismember(k,peaks_glove(:,1))     %for thumb            
                thumb=1;
            end
            if ismember(k,peaks_glove(:,2))      %for index finger
                index=1; 
            end
            
            if ismember(k,peaks_glove(:,4))     %for ring finger
                ring=1;
            end
            if ismember(k,peaks_glove(:,5))     %for little finger
                little=1;
            end
            k=k+1;
        end
        if   middle==1  && thumb==0 && index==0 && ring ==0&& little ==0  %for middle
            data_num=data_num+1;
            disp(["Time Segment No (Middle):" num2str(data_num)]);
            pause(2);
            a=start-2000;
            middle_data(:,:,data_num)=ecog_data(a:a+3999,:);           
        end
        i=i+2099 ;
    end
    i=i+1;
    
end

%% Temporal feature extraction for middle
temporal_middle=[];

n=size(middle_data,3);
for i=1:n
    disp(["Feature Extraction for Time Segment No (Middle):" num2str(i)]);
    pause(2);
    channel=[];
    for j=1:num_channels                %for each channel
        channel_task=[];
        channel_relax=[];
        
        start=1;
        for k=1:76                    %take periodogram for relax time
            window=middle_data(start:start+249,j,i); %250ms window    %%%%%%only change here
            [Pxx_relax,F_relax] = periodogram(window,hamming(length(window)),length(window),Fs);
            start=start+50;   %stride 50 ms
            Pxx_relax=10*log10(Pxx_relax);
            channel(j,k)=mean(Pxx_relax(low_cut_1:upper_cut_1));        %(2:18) 2:3 3:4 4:11 11:18 19:34 35:51
            channel(j+num_channels,k)=mean(Pxx_relax(low_cut_2:upper_cut_2));   %(19:76)
            
        end
        
    end
    
    temporal_middle(:,:,i)=channel;
    
end


%% data of Middle

channel_mean=mean(temporal_middle(:,1:16,:),2);   %1:16 before
temporal_middle_normalized=temporal_middle-channel_mean;
temporal_middle_normalized_filt= imfilter(temporal_middle_normalized,fspecial('gaussian',[3,3],1));
temporal_middle_normalized_all=mean(temporal_middle_normalized_filt,3);



figure;
f=1:num_channels;
t=0.125:0.05:3.875;
imagesc(t,f,temporal_middle(1:num_channels,:,1));   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("before nomalizing Middle Finger")
figure;
imagesc(t,f,temporal_middle_normalized_all(1:num_channels,:));
grid on
colormap('jet')
title("after nomalizing Middle Finger 4-8 Hz data name : wm")
figure;
imagesc(t,f,temporal_middle_normalized_all(num_channels+1:num_channels*2,:));
grid on
colormap('jet')
title("after nomalizing Middle Finger 8-12 Hz data name : wm")

%% for Ring

ecog_data=notch_filtered();
relax_times=extracttask(stim,4);   %%% Ring==4
sum(relax_times)/2000
[n,m]=size(ecog_data);
ring_data=[];                     %change here
data_num=0;
i=1;
while i<=n-1900
    if relax_times(i)~=0
        k=i+299;
        first=0;
        thumb =0;
        index=0;
        middle=0;
        ring=0;
        little=0;
        while k<=i+1699                         % checking whther the subject has done the task correctly or not using glove signal
            if ismember(k,peaks_glove(:,4))     %for ring finger
                ring=1;
                if first==0
                    start=k;
                    first=1;
                end
                break
            end
            if ismember(k,peaks_glove(:,1))     %for thumb            
                thumb=1;
            end
            if ismember(k,peaks_glove(:,2))      %for index finger
                index=1; 
            end
            if ismember(k,peaks_glove(:,3))     %for middle finger
                middle=1;
            end
            
            if ismember(k,peaks_glove(:,5))     %for little finger
                little=1;
            end
            k=k+1;
        end
        if   ring ==1 && thumb==0 && index==0 && middle==0 %&& little ==0   %for ring
            data_num=data_num+1;
            disp(["Time Segment No (Ring):" num2str(data_num)]);
            pause(2);
            a=start-2000;
            ring_data(:,:,data_num)=ecog_data(a:a+3999,:);           
        end
        i=i+2099 ;
    end
    i=i+1;
    
end

%% Temporal feature extraction for ring
temporal_ring=[];

n=size(ring_data,3);
for i=1:n
    disp(["Feature Extraction for Time Segment No (Ring):" num2str(i)]);
    pause(2);
    channel=[];
    for j=1:num_channels                %for each channel
        channel_task=[];
        channel_relax=[];
        
        start=1;
        for k=1:76                    %take periodogram for relax time
            window=ring_data(start:start+249,j,i); %250ms window    %%%%%%only change here
            [Pxx_relax,F_relax] = periodogram(window,hamming(length(window)),length(window),Fs);
            start=start+50;   %stride 50 ms
            Pxx_relax=10*log10(Pxx_relax);
            channel(j,k)=mean(Pxx_relax(low_cut_1:upper_cut_1));        %(2:18) 2:3 3:4 4:11 11:18 19:34 35:51
            channel(j+num_channels,k)=mean(Pxx_relax(low_cut_2:upper_cut_2));   %(19:76)
            
        end
        
    end
    
    temporal_ring(:,:,i)=channel;
    
end


%% data of Ring

channel_mean=mean(temporal_ring(:,1:16,:),2);
temporal_ring_normalized=temporal_ring-channel_mean;
temporal_ring_normalized_filt= imfilter(temporal_ring_normalized,fspecial('gaussian',[3,3],1));
temporal_ring_normalized_all=mean(temporal_ring_normalized_filt,3);



figure;
f=1:num_channels;
t=0.125:0.05:3.875;
imagesc(t,f,temporal_ring(1:num_channels,:,1));   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("before nomalizing Ring Finger")
figure;
imagesc(t,f,temporal_ring_normalized_all(1:num_channels,:));
grid on
colormap('jet')
title("after nomalizing Ring Finger 4-8 Hz data name : wm")
figure;
imagesc(t,f,temporal_ring_normalized_all(num_channels+1:num_channels*2,:));
grid on
colormap('jet')
title("after nomalizing Ring Finger 8-12 Hz data name : wm")


%% for Little

ecog_data=notch_filtered();
relax_times=extracttask(stim,5);   %%% little==4
sum(relax_times)/2000
[n,m]=size(ecog_data);
little_data=[];                     %change here
data_num=0;
i=1;
while i<=n-1900
    if relax_times(i)~=0
        k=i+299;
        first=0;
        thumb =0;
        index=0;
        middle=0;
        ring=0;
        little=0;
        while k<=i+1699                         % checking whther the subject has done the task correctly or not using glove signal
            if ismember(k,peaks_glove(:,5))     %for little finger
                little=1;
                if first==0
                    start=k;
                    first=1;
                end
                break
            end
            if ismember(k,peaks_glove(:,1))     %for thumb            
                thumb=1;
            end
            if ismember(k,peaks_glove(:,2))      %for index finger
                index=1; 
            end
            if ismember(k,peaks_glove(:,3))     %for middle finger
                middle=1;
            end
            if ismember(k,peaks_glove(:,4))     %for ring finger
                ring=1;
            end
            
            k=k+1;
        end
        if   little ==1  && thumb==0 && index==0 && middle==0% && ring ==0  %for little
            data_num=data_num+1;
            disp(["Time Segment No (Little):" num2str(data_num)]);
            pause(2);
            a=start-2000;
            little_data(:,:,data_num)=ecog_data(a:a+3999,:);           
        end
        i=i+2099 ;
    end
    i=i+1;
    
end

%% Temporal feature extraction for little
temporal_little=[];

n=size(little_data,3);
for i=1:n
    disp(["Feature Extraction for Time Segment No (Little):" num2str(i)]);
    pause(2);
    channel=[];
    for j=1:num_channels                %for each channel
        channel_task=[];
        channel_relax=[];
        
        start=1;
        for k=1:76                    %take periodogram for relax time
            window=little_data(start:start+249,j,i); %250ms window    %%%%%%only change here
            [Pxx_relax,F_relax] = periodogram(window,hamming(length(window)),length(window),Fs);
            start=start+50;   %stride 50 ms
            Pxx_relax=10*log10(Pxx_relax);
            channel(j,k)=mean(Pxx_relax(low_cut_1:upper_cut_1));        %(2:18) 2:3 3:4 4:11 11:18 19:34 35:51
            channel(j+num_channels,k)=mean(Pxx_relax(low_cut_2:upper_cut_2));   %(19:76)
            
        end
        
    end
    
    temporal_little(:,:,i)=channel;
    
end


%% data of little

channel_mean=mean(temporal_little(:,1:16,:),2);
temporal_little_normalized=temporal_little-channel_mean;
% 
% channel_filt= imfilter(channel,fspecial('gaussian',[20 20],1));
%temporal_little_normalized_filt= imfilter(temporal_little_normalized,fspecial('gaussian',[5,5],1));
temporal_little_normalized_filt= imfilter(temporal_little_normalized,fspecial('gaussian',[3,3],1));
temporal_little_normalized_all=mean(temporal_little_normalized_filt,3);



figure;
f=1:num_channels;
t=0.125:0.05:3.875;
imagesc(t,f,temporal_little(1:num_channels,:,1));   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("before nomalizing Little Finger")
figure;
imagesc(t,f,temporal_little_normalized_all(1:num_channels,:));
grid on
colormap('jet')
title("after nomalizing Little Finger  4-8 Hz data name : wm")
figure;
imagesc(t,f,temporal_little_normalized_all(num_channels+1:num_channels*2,:));
grid on
colormap('jet')
title("after nomalizing Little Finger 8-12 Hz data name : wm")


%% Save Extracted Features
save(strcat("./extracted_features/Thumb_temporal_", sub_name, "_", freq_band_1, "Hz_", freq_band_2, "Hz_ham"), "temporal_thumb_normalized");
save(strcat("./extracted_features/Index_temporal_", sub_name, "_", freq_band_1, "Hz_", freq_band_2, "Hz_ham"), "temporal_index_normalized");
save(strcat("./extracted_features/Middle_temporal_", sub_name, "_", freq_band_1, "Hz_", freq_band_2, "Hz_ham"), "temporal_middle_normalized");
save(strcat("./extracted_features/Ring_temporal_", sub_name, "_", freq_band_1, "Hz_", freq_band_2, "Hz_ham"), "temporal_ring_normalized")
save(strcat("./extracted_features/Little_temporal_", sub_name, "_", freq_band_1, "Hz_", freq_band_2, "Hz_ham"), "temporal_little_normalized")
