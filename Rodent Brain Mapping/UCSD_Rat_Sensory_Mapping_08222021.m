% Sensory mapping codes used for rat whisker barrel mapping
% Last updated by Youngbin Tchoe 08/23/2021
%% (1) Read the impedance information
% Load the electrode impedance data (.csv)
imp_filepath = '.\Rat2_Impedance_in_saline.csv'; % input the file location
% Channel selection criteria: 1 kHz impedance magnitude
lower_bound = 1e3; % larger than 1 kOhm
upper_bound = 100e3; % lower than 100 kOhm
% Call the impedance mapping function
k = orb1024_impedance_ratHD_yb(imp_filepath, lower_bound, upper_bound);
    % k vector will save channels with good impedance (> 1 kOhm && < 100 kOhm)

%% (2) Read Intan .rhd file
% Load the .rhd recording data
rhd_filepath = '.\Left_Whisker_E4_blue_200317_145821.rhd'; % input the file location
% Call the function to unfold the data - This is a customized one
% The function provided by Intan can be found at: 
% https://intantech.com/downloads.html?tabSelect=Software&yPos=0
read_Intan_RHD2000_file_custom(rhd_filepath)
    % time: t_amplifier
    % potential of 1024 channels: amplifier_data
    % analog input data: boad_adc_data
% Save sampling rate
fs = frequency_parameters.amplifier_sample_rate;

%% (3) Load electrode mapping
    load('elecCoords_ratHD');
    elecX = elecCoords_ratHD(:,1);
    elecY = elecCoords_ratHD(:,2);
    
%% Generate 4th order Butterworth filter
    N = 4; % 4th order butter worth
    lf = 70; % lower bound in Hz
    hf = 190; % Upper bound in Hz
    % Define butterworth filter
    d_HG = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf,'HalfPowerFrequency2',hf, ...
    'SampleRate',fs);

%% Capturing the trials

% Setup the time segment
timeStart = -0.2; % (second) How much time before the stimulus
timeEnd = 0.8; % (second) How much time after the stimulus
timeDur = -timeStart + timeEnd; % (second) Time duration captured for each trial

% N: Total data length
N = length(board_adc_data(1,:));

% TTL time stamp capture from the Analog input
findPeak = 1; % Loop captures trial only when this is 1
m_temp = 0; % Saves the last captured index number
trialIndex = []; % Index numbers of each trial are saved here
adc_time = 0.5; % (second) % This time should be longer than a single pulse captured in Analog input
for m = 1:N
    
    if(m - m_temp) > adc_time*fs; % make findPeak = 1 after a single pulse passed
        findPeak = 1;
    end
    
    % If the amplitude is above 2.5V, save the index number
    if ( board_adc_data(1,m) > 2.50) && (findPeak)
       findPeak = 0; 
       m_temp = m;
       trialIndex = [trialIndex m];
    end

end

% Removing the last trial if it is not fully recorded
if (t_amplifier(end) - trialIndex(end)/fs) < (timeEnd)
	numTrials = length(trialIndex)-1;
else
    numTrials = length(trialIndex);
end

% Display number of trials
display(sprintf("Number of trial: %d", numTrials));


%% Trial averaging using captured index numbers

ampDataTA = zeros(1024,fs+1); % Trial averaged data will be stored here
for i = 1:numTrials
index1 = trialIndex(i) + timeStart*fs; % Index number 1, -0.2 s before stim
index2 = trialIndex(i) + timeEnd*fs; % Index number 2, +0.8 s after stim
ampDataTA = ampDataTA + amplifier_data(:, index1:index2); % Summation
display(sprintf("Summation of trial number %d out of %d trials", i, numTrials));
end
ampDataTA = ampDataTA/numTrials; % Devide by the number of trial
display(sprintf("Trial averaging done"));

%% Obtain the RMS of the response at a specific time range: 5 ms ~ 100 ms after stim
ampDataTA_RMS = zeros(1024,2); 
time_ms1 = 5; % 5 ms after stim
time_ms2 = 100; % 50 ms after stim
range_ms = int32((abs(timeStart)+time_ms1*1e-3)*fs):int32((abs(timeStart)+time_ms2*1e-3)*fs);

for i = 1:1024
	ampDataTA_RMS(i,1) = rms(ampDataTA(i,range_ms)); % RMS data
    ampDataTA_RMS(i,2) = i; % Channel number
end
display(sprintf("Obtained RMS data"));

%% Common mode subtraction using channels with low responses

    % Sorting the channels according to the RMS value
    % 'Y_min' is the vector sorted according to the RMS value
    X = [ampDataTA_RMS(k,1) ampDataTA_RMS(k,2)];
    Y_min = sortrows(X, 'ascend');

    % Taking the common mode
    num_CM = 10; % Number of minimum RMS channels to take common mode
    ampDataCM = zeros(1024,fs+1); % Common mode data will be stored here
    
    % Averaging the 10 lowest RMS channels
    for i = 1:num_CM
    ampDataCM = ampDataCM + ampDataTA(Y_min(i, 2),:); 
    end
    ampDataCM = ampDataCM/num_CM;

    % Common-mode rejected data: ampDataTaCMR
    ampDataTaCMR = zeros(1024,fs+1);
    ampDataTaCMR = ampDataTA - ampDataCM;
    
    display(sprintf("Obtained Common-mode rejected data"));
    
%% RMS of Trial averaged & Common mode rejection signals
    ampDataTaCMR_RMS = zeros(1024,2);
    for i = 1:1024
            ampDataTaCMR_RMS(i,1) = rms(ampDataTA(i,:));
            ampDataTaCMR_RMS(i,2) = i;
    end
    X = [ampDataTaCMR_RMS(k,1) ampDataTaCMR_RMS(k,2)];
    Y = sortrows(X, 'ascend');
    
    display(sprintf("Obtained Common-mode rejected RMS data"));
    
    % Plot the common mode rejected - RMS mapping data
    size = 500; % Size of the figure
%     [M, posMaxIndex] = max(ampDataTaCMR_RMS(k,1)) % Finding max RMS
    [m, posMinIndex] = min(ampDataTaCMR_RMS(k,1)) % Finding min RMS
    
    % Plot range
    v_max = 5*std(ampDataTaCMR_RMS(k,1)); % 5 sigma : Maximum voltage (uV)
    v_min = m; % uV - min
    figure1=figure('Position', [150, 150, 150+size, 150+size]);
    
    % Scatter plot
    scatter(elecX(k),elecY(k),200,ampDataTaCMR_RMS(k,1),'filled', 's'); 
    
    % Plot parameters
    box on;
    caxis([v_min v_max]); % Color range
    colorbar % Enable colorbar
    colormap('jet'); % Color map
    pbaspect([32 32 1]); % Aspect ratio of the plot
    axis([-2.5 2.5 -2.5 2.5]); % Plot range
    xlabel('mm'); % x label
    ylabel('mm'); % y label
    txtTitle = append('RMS_mapping');
    title(txtTitle,'interpreter', 'none');
    h = colorbar;
    ylabel(h, 'Voltage (\muV)');
    set(gca,'FontSize',12);
    
	display(sprintf("RMS mapping was plotted"));
    
%% High gamma activity RMS mapping
prefix = 'HGA';

t0 = cputime; % Store time
data_HG = [];
data_HG_Amp = [];

for i = 1:1024

    % 60Hz and harmonics Bandstop
    data_60 = ampDataTaCMR(i, :);

     for nf = [60,120,180,240,300,360,420,480]
        Wo = nf/(fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        data_60 = filtfilt(b,a,data_60);
     end
     
% High gamma bandpass filtering 70~190 Hz 
data_HG(i,:) = filtfilt(d_HG,data_60);
 
% Amplitude of the data: Hilbert transform & abs
data_HG_Amp(i,:) = abs(hilbert(data_HG(i,:)));

t1 = cputime;
display(sprintf("High gamma filtering: %u percent \t %d out of 1024 completed \t %.0f sec elapsed", uint8(100*i/1024), i, (t1-t0)));
end

display(sprintf("High gamma filetering is completed"));

%% High gamma RMS plot

% Select the time range to take RMS 
time1 = 0.020 % (second) from 20 ms 
time2 = 0.050 % (secod) to 50 ms
range = int32((time1+0.2)*fs+1):int32((time2+0.2)*fs+1);

% Obtain RMS data of high gamma
data_HG_RMS = zeros(1024,1);
for i = 1:1024
        data_HG_RMS(i) = rms(data_HG(i,range));
end

% Finding peak RMS location
% [M, posMaxIndex_HG] = max(data_HG_RMS(k))
% [m, posMinIndex_HG] = min(data_HG_RMS(k))
% Plot the RMS mapping data
size = 500; % size of the plot
v_max = 5*std(data_HG_RMS(k)); % 5 sigma
v_min = 0;

figure2=figure('Position', [150, 150, 150+size, 150+size]);

% Scatter mapping of high gamma activity
scatter(elecX(k),elecY(k),180,data_HG_RMS(k),'filled', 's'); 

% Plot parameters
box on;
caxis([v_min v_max]); 
colorbar
colormap(jet);
pbaspect([32 32 1]);
axis([-2.5 2.5 -2.5 2.5]);
xlabel('mm');
ylabel('mm');
txtTitle = append('HGA_RMS_', num2str(int32(time1*1E3)), 'to', num2str(int32(time2*1E3)),'ms');
title(txtTitle,'interpreter', 'none');
h = colorbar;
ylabel(h, 'Voltage (\muV)');
set(gca,'FontSize',12);

    %% Mapping plot of waveforms
    % Downsampling to 1 kSamples/s
	fs_ds = 1e3; % Hz
    raw_TA_data = downsample(ampDataTaCMR', fs/fs_ds)';
    t_ms1 = 0; % from 0 ms after stim
    t_ms2 = 80; % to 80 ms after stim
    range = (-timeStart+t_ms1*1e-3)*fs_ds:(-timeStart+t_ms2*1e-3)*fs_ds;
    size = 1000; % Figure size
    figure3=figure('Position', [0, 0, size, size]);
    
    % Coloration of wavetraces based on their RMS power
    raw_TA_rms = [];
    for i = 1:1024
        raw_TA_rms(i) = rms(raw_TA_data(i,range));
    end
	raw_TA_rms = raw_TA_rms/max(raw_TA_rms(k));

    % Plotting the wavetraces
    x_offset = -0.075;
    time_adj = 600; % Should find appropriate contant to set wavetrace length
    volt_adj = 2000; % Should find appropriate contant to set wavetrace height
    
    for i = k
        plot(elecX(i)+x_offset+(1:length(range))/time_adj,elecY(i)+raw_TA_data(i,range)/volt_adj, 'Color',[raw_TA_rms(i),0,0], 'LineWidth', 1);
        hold on;
    end
    
    % Plot parameters
    xlabel('Position (mm)');
    ylabel('Position (mm)');
    pbaspect([32 32 1]);
    axis([-2.5 2.5 -2.5 2.5]);
    line_width = 1.5;
    font_size = 15;
    set(gca,'FontSize',font_size);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',line_width);

%% Waterfall plot

    t_ms1 = -50; % from -50 ms before stim
    t_ms2 = 150; % to 150 ms after stim
    range = int32((-timeStart+t_ms1*1e-3)*fs_ds):int32((-timeStart+t_ms2*1e-3)*fs_ds);
    size = 1000; % plot size

    figure4=figure('Position', [0, 0, size/2, size]);

    dist_sort = zeros(1024,2); % Sort the channel according to the distance from the (cent_x, cent_y)
    cent_x = -2; % Center position, x (mm)
    cent_y = -2; % Center position, y (mm)
    v_offset = 2; % Voltage offset between the traces
    
    for i = 1:1024 % Sorting distance
            dist_sort(i,1) = sqrt((elecX(i)-cent_x)^2+(elecY(i)-cent_y)^2);
            dist_sort(i,2) = i;
    end
    X = [dist_sort(k,1) dist_sort(k,2)];
    Y = sortrows(X, 'ascend');

    for i = 1:2:length(Y)-1 % Plot all traces
        plot((1:length(range))-50, i*v_offset+raw_TA_data(Y(i,2),range));
        plot((1:length(range))-50,-i*v_offset+raw_TA_data(Y(i+1,2),range));
        hold on;
    end
    
    % Plot parameters
    pbaspect([16 32 1]);
    xline(0);
    ylim([-v_offset*(length(Y)+10)-20 v_offset*(length(Y)+10)]);
    xlim([-50 100]);
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');

