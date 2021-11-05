% read_Intan_RHD2000_file
load('elecCoords_8x8_R.mat');
% save('CS_BPF.mat', 'SOS_CS', 'G_CS');
load('CS_BPF.mat');
 elecX = elecCoords_8x8_R(:,1);
 elecY = elecCoords_8x8_R(:,2);
% load('k_vector_Sep08_2020.mat');
k = k_R;
ch_num = 1024;
% save('k_vector_Sep04_2020.mat', 'k');

% plot(board_adc_data)
% plot(t_amplifier)
% plot(t_amplifier, board_adc_data)
% 
% i = 1;
% 
% for i = 1:100:1024
% plot(t_amplifier,5*i+amplifier_data(i,:))
% hold on;
% end
%%
     % Filtering
    fs = 2e4;
    N = 4;
    lf = [1 4 8 13 30 70 30];% 190 250 500 750]; % lower bound in Hz
    hf = [4 8 12 30 70 190 3000];% 250 500 750 2000]; % Upper bound in Hz
    
    i_freq = 7;% 7: Central Sulcus
    d_CS = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
    'SampleRate',fs);

%% Stimulation time stamp capture
% detection by amplitude
isPeak = 1;
k_temp = 0;
pkIndex = [];
fs = 20000;
t1 = 560%190; % stim start (s)
t2 = 585%201; % stim end (s)
th = 50;
index1 = int32(fs*(t1 - t_amplifier(1)));
index2 = int32(fs*(t2 - t_amplifier(1)));
% index1 = 1.05 * 1e6;
% index2 = 1.20 * 1e6;

x = highpass(amplifier_data(30,:), 1000, fs);
plot(t_amplifier, x)

% plot(t_amplifier(:),x)
for m = index1:index2
    % Minimum time distance between TTLs
    if(m - k_temp) > 0.1*fs;
        isPeak = 1;
    end
    % If the amplitude is above 0.1V, save
    if ( x(1,m) > th) && (isPeak)
       isPeak = 0; 
       k_temp = m;
       pkIndex = [pkIndex m];
    end
end
length(pkIndex)
% 
figure;
plot(t_amplifier, x)

hold on;
plot(t_amplifier(pkIndex), th*ones(length(pkIndex)), '-o')

plot(x)
hold on;
plot(pkIndex, th*ones(length(pkIndex)), '-o')


%% 'k' saves the channels in reliable impedance magnitude range and excludes
% k = [];
% j = 1;
%  for i = 1:1024
%      if(impedanceMag(i)<100e3 && impedanceMag(i)>1e3) %&& (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1))&& (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1))&& (~(abs(elecX(i)+1.425) < 0.1 && abs(elecY(i)+2.325) < 0.1)) && (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1)) && (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1)) && (~((abs(elecX(i)-0.15) < 0.15) && (elecY(i)>2.3))) && (~(abs(elecX(i)-0.525) < 0.1 && abs(elecY(i)-2.325) < 0.1)) && (~(abs(elecX(i)-2.325) < 0.1 && abs(elecY(i)-2.175) < 0.1)) && (~(abs(elecX(i)-2.175) < 0.1 && abs(elecY(i)-2.175) < 0.1))
%          k(1,j) = i;
%          j = j+1;
%      end
%  end
%  length(k)
%  save('k_vector_Sep04_2020.mat', 'k');
% load('k_vector_Sep04_2020.mat');
%% Averaging one second time segment data into an array
timeStart = -0.1;
timeInterval = 0.5;
numTrials = length(pkIndex);

% Real trial averaging loop
ampDataTA = zeros(1024,fs*timeInterval+1);
for i = 1:numTrials-1
index1 = pkIndex(i) + timeStart*fs; % time1
index2 = pkIndex(i) + (timeStart+timeInterval)*fs; % time2

    y_bpf = zeros(1024, length(index1:index2));
    for j = 1:1024
    x = amplifier_data(j, index1:index2);
    y_bpf(j,:) = filtfilt(SOS_CS,G_CS,x);
    end

ampDataTA = ampDataTA + y_bpf;

% plot((1:length(index1:index2))/fs, amplifier_data(k(:), index1:index2));
i
end
ampDataTA = ampDataTA/numTrials;

% plot((1:length(index1:index2))/fs, ampDataTA(k(:), :));

%     % 60Hz and harmonics Bandstop
%     data_60 = ampDataTaCMR(i, :);
%       
%      for nf = [60,120,180,240,300,360,420,480]
%         Wo = nf/(fs/2);
%         BW = Wo/35;
%         [b,a] = iirnotch(Wo, BW);
%         data_60 = filtfilt(b,a,data_60);
%      end

% 
% for i = 1:100:1024
% plot((1:length(ampDataTA))/fs,0.1*i+ampDataTA(i,:))
% hold on;
% end

    %% Sort channels
    offset = 10;
    fs = 2e4;
    

    for x_pos = 0.9:1.8:40.5%-40.5:1.8:0%-15.5:1:15.5;%0.9:1.8:40.5
        
	k_temp = [];
    m = 1;
    for i = 1:length(k)
        if abs(elecX(k(i)) - x_pos) < 0.1
            k_temp(m,1) = elecY(k(i));
            k_temp(m,2) = k(i);
            m = m+1;
        end
    end
    
    if isempty(k_temp)
        continue;
    end
    
    X = [k_temp(:,1) k_temp(:,2)];
    k_sorted = sortrows(X, 'ascend');
    
    pos_xy = 10;
    window_size_x = 500;
    window_size_y = 1000;
    figure1=figure('Position', [pos_xy, pos_xy, pos_xy+window_size_x, pos_xy+window_size_y]);
    for i = 1:length(k_temp)
    plot((1:length(ampDataTA))/fs+timeStart,offset*i+ampDataTA(k_sorted(i,2),:))
    hold on;
    end
    xlim([-0.05 0.1]);
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ylim([-offset offset*(length(k_temp)+2)]);
    txt = ['Time = ' num2str(x_pos, '%.1f mm')];
    title(txt);
    
    
    txt_filename = append('CS_2048L_trim_LUstim_Sep8_2020_', num2str(t1), 's_', num2str(x_pos), 'mm_', '.png');;
    print('-dpng','-r600',txt_filename);
    
    close(figure1);
    
    end
    
   
 %% Mapping plot of ampDataTaCMR
 
 % save('SSEP_2048L_LU.mat', 'ampDataTA', 'k_L');
 % save('SSEP_2048R_LU.mat', 'ampDataTA', 'k_R');
    load('SSEP_2048R_LU.mat');
    k = k_R;
    load('elecCoords_8x8_R.mat');
    elecX = elecCoords_8x8_R(:,1);
    elecY = elecCoords_8x8_R(:,2);
    
    fs = 2e4;
    raw_TA_data = downsample(ampDataTA', fs/1000)';
    fs = 1e3;
    range = 0.105*fs:0.150*fs;
	x_size = 1200;
    y_size = 1200;
    line_width = 1;
    font_size = 10;
    figure1=figure('Position', [0, 0, x_size, y_size]);
    t_div = 30;
    v_div = 15;
    
    raw_TA_rms = [];
    for i = 1:1024
        raw_TA_rms(i) = rms(raw_TA_data(i,range));
    end
    plot(raw_TA_rms(k))
    raw_TA_rms = raw_TA_rms/max(raw_TA_rms(k));
%     scatter(elecX(k),elecY(k),150,raw_TA_data(k,142),'filled', 's'); 
%     colorbar;
    
%     cm = colormap(jet); 
    t_peak_ds = 100+23; %ms
    l_width = 0.5;
    max_amp = 16; %max(max(abs(raw_TA_data(k,t_peak_ds))));
    
    for i = k
        if raw_TA_data(i,t_peak_ds) > 0
        plot(elecX(i)-0.075+(1:length(range))/t_div, elecY(i)+raw_TA_data(i,range)/v_div, 'Color',[raw_TA_data(i,t_peak_ds)/max_amp,0,0], 'LineWidth', l_width);
        else
        plot(elecX(i)-0.075+(1:length(range))/t_div, elecY(i)+raw_TA_data(i,range)/v_div, 'Color',[0,0,-raw_TA_data(i,t_peak_ds)/max_amp], 'LineWidth', l_width);
        end
        hold on;
    end
    
    xlabel('Position (mm)');
    ylabel('Position (mm)');
    pbaspect([x_size y_size 1]);
    lim = 17   ;% - 0.075;
    axis([0 42 -42 42]);
%     axis off;

    set(gca,'FontSize',font_size);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',line_width);
    
    txt = 'CS_mapping_2048R_Sep08_2020';
    print('-dpng','-r300',txt);