% read_Intan_RHD2000_file
load('elecCoords_CS_32x32');
load('k_vector_Linear_CS_11132020.mat');
% save('CS_BPF.mat', 'SOS_CS', 'G_CS');
load('CS_BPF.mat');
 elecX = elecCoords_CS_32x32(:,1);
 elecY = elecCoords_CS_32x32(:,2);
 
     % Filtering
    fs = 2e4;
    N = 4;
    lf = [1 4 8 13 30 70 30];% 190 250 500 750]; % lower bound in Hz
    hf = [4 8 12 30 70 190 3000];% 250 500 750 2000]; % Upper bound in Hz
    
    i_freq = 7;% 7: Central Sulcus
    d_CS = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
    'SampleRate',fs);

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

%% Stimulation time stamp capture
% detection by amplitude
isPeak = 1;
k_temp = 0;
pkIndex = [];
fs = 20000;
t1 = 130%157%140%66%181; % stim start (s)
t2 = 140%170%170%99%210; % stim end (s)
th = 15;
index1 = fs*int32(t1 - t_amplifier(1));
index2 = fs*int32(t2 - t_amplifier(1));

x = highpass(amplifier_data(32,:), 5000, fs);
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

%% 'k' saves the channels in reliable impedance magnitude range and excludes
% k = [];
% j = 1;
%  for i = 1:1024
%      
%      if(impedanceMag(i)<100e3 && impedanceMag(i)>1e3) 
%          %&& (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1))&& (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1))&& (~(abs(elecX(i)+1.425) < 0.1 && abs(elecY(i)+2.325) < 0.1)) && (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1)) && (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1)) && (~((abs(elecX(i)-0.15) < 0.15) && (elecY(i)>2.3))) && (~(abs(elecX(i)-0.525) < 0.1 && abs(elecY(i)-2.325) < 0.1)) && (~(abs(elecX(i)-2.325) < 0.1 && abs(elecY(i)-2.175) < 0.1)) && (~(abs(elecX(i)-2.175) < 0.1 && abs(elecY(i)-2.175) < 0.1))
%          k(1,j) = i;
%          j = j+1;
%      end
%  end
%  length(k)
%  save('k_vector_Linear_CS_11132020.mat', 'k');
load('k_vector_Linear_CS_11132020.mat');
%% Averaging one second time segment data into an array
timeStart = -0.1;
timeInterval = 0.5;
numTrials = length(pkIndex);

% Real trial averaging loop
ampDataTA = zeros(1024,fs*timeInterval+1);
    y_bpf = zeros(numTrials, 1024, int32(timeInterval*fs)+1);

for i = 1:numTrials
index1 = pkIndex(i) + timeStart*fs; % time1
index2 = pkIndex(i) + (timeStart+timeInterval)*fs; % time2

    for j = 1:1024
    x = amplifier_data(j, index1:index2);
%       y_bpf(j,:) = filtfilt(SOS_CS,G_CS,x);
        y_bpf(i, j,:) = reshape(filtfilt(d_CS,x), [1,1,10001]);

    end
ampDataTA = ampDataTA + reshape(y_bpf(i,:,:), [1024, 10001]);
% ampDataTA = ampDataTA + amplifier_data(:, index1:index2);%y_bpf;

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

%% Referencing: Substracting minimum RMS 
% % Common mode subtraction from minimum RMS channel
% % RMS of the response
% 
%         t1 = 27;
%     timeStart = -5;
%     ampDataTA = amplifier_data(:, (t1+timeStart)*fs:(t1+timeStart+2*abs(timeStart))*fs);
% 
%     ampDataTA_RMS = [];
%     for i = 1:1024
%             ampDataTA_RMS(i,1) = rms(ampDataTA(i,:));
%             ampDataTA_RMS(i,2) = i;
%             i
%     end
%     X = [ampDataTA_RMS(k,1) ampDataTA_RMS(k,2)];
%     Y = sortrows(X, 'ascend');
%     
%     num_CM = 10;%length(k); % Number of minimum RMS channels to take common mode
%     timeInterval = 10;
%     ampDataCM = zeros(1024,fs*timeInterval+1);
% 
%     for i = 1:num_CM
%     ampDataCM = ampDataCM + ampDataTA(Y(i, 2),:);
%     i*1e-3
%     end
%     ampDataCM = ampDataCM/num_CM;
% 
%     ampDataTaCMR = zeros(1024,fs*timeInterval+1);
%     ampDataTaCMR = ampDataTA - ampDataCM;
% %     
% %     for i = 1:100:1024
% %     plot((1:10001)/fs,0.1*i+ampDataTaCMR(i,:))
% %     hold on;
% %     end
   
    %% Sort channels
    fs = 2e4;
    data = ampDataTA;

    for x_pos = 2.5%-15.5:1:15.5;
        
	k_temp = [];
    m = 1;
    for i = 1:length(k)
        if abs(elecX(k(i)) - x_pos) < 0.1
            k_temp(m,1) = elecY(k(i));
            k_temp(m,2) = k(i);
            m = m+1;
        end
    end
    X = [k_temp(:,1) k_temp(:,2)];
    k_sorted = sortrows(X, 'ascend');
    
    pos_xy = 10;
    window_size_x = 500;
    window_size_y = 1000;
    figure1=figure('Position', [pos_xy+500, pos_xy-50, pos_xy+window_size_x, pos_xy+window_size_y]);
   
%     data = smoothdata(ampDataTA, 'gaussian', 50);
% 	data = ampDataTA;
% 	data = ddf;
	offset = 10*std(reshape(data, [1024*length(data),1]));
    num_ln = length(k_temp);
      for i = 1:num_ln
        t = 1e3*((1:length(data(1,:)))/fs+timeStart);
%                 yline(offset*i, 'b');
%                 hold on;
        plot(t,offset*i+data(k_sorted(i,2),:), 'k', 'Linewidth',1.25)   
        hold on;
      end
%       set(gcf,'GraphicsSmoothing','off') 
      
%     i = 9;
%     range = int32((0.1+0.001)*fs):int32((0.1+0.006)*fs);
% 	plot((1:length(data(k_sorted(i,2),range)))/fs-0.1, data(k_sorted(i,2),range), 'k')   
%     rms(data(k_sorted(i,2),range))
    
%     70/1.0895

    
%     for m = 1:numTrials
% 	data = reshape(y_bpf(m,:,:), [1024, 10001]);
%     for i = 1:length(k_temp)
%         t = 1e3*((1:length(data(1,:)))/fs+timeStart);
%         plot(t,offset*i+data(k_sorted(i,2),:), 'color', [m/numTrials,0,0])      
%         hold on;
%     end
%     end
    box on;
% 	xline(20, 'r');
% 	xline(25, 'r');
%     xline(28, 'r');
% 	xline(35);
%     xlim([-0.05 0.1]);
	xlim([-50 100]);

    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    ylim([-offset offset*(num_ln+2)]);
    txt = ['X = ' num2str(x_pos, '%.1f mm')];
    title(txt);
    
    
    txt_filename = append('CS3_Nov13_2020_', num2str(t1), 's_', num2str(x_pos), 'mm_', '.png');;
    print('-dpng','-r600',txt_filename);
    
%     close(figure1);
    
    end
    
    %%
    window_size_x = 800;
    window_size_y = 800;
    figure1=figure('Position', [pos_xy+500, pos_xy-50, pos_xy+window_size_x, pos_xy+window_size_y]);
   
    t_point = int32((0.1+0.028)*fs);
	plot(data(k_sorted(1:32,2),t_point), '-o', 'linewidth', 1.5);
    yl = yline(0);
	set(yl,'linewidth', 1.5);

%     xline(15);
    xlim([ 0 33]);
    xlabel('Channel');
    ylabel('Potential (\muV)');
    ax=gca;ax.LineWidth=1.5;
    set(gca,'fontsize', 15);
        set(gca,'fontweight', 'bold');
    txt_filename = 'SSEP_along_x0p5.png';
%     print('-dpng','-r300',txt_filename);
    %% 
    df = [];
    ddf = [];
    for i = 1:1024
       
        f = smoothdata(ampDataTA(i,:), 'gaussian', 50);
        df(i,:) = diff(f);
        ddf(i,:) = diff(diff(f));
        i
        
    end
        %% Data
    data = [];
    for i = 1:1024
        data(i,:) = smoothdata(ampDataTA(i,:), 'gaussian', 50);
        i
    end
    %% Plot at specific time
%         data = ampDataTA;
%         data = smoothdata(ampDataTA, 'gaussian', 50);
%         data = ddf;
        
    for time_ms = 28%[10:1:50, 60:10:100]
    ref_ms = 25;
    time = time_ms * 1e-3;
    specific_index = int32((time-timeStart)*fs);
    ref_index = int32((ref_ms*1e-3-timeStart)*fs);
%         REF = zeros(length(k),1);
        REF = zeros(1024,1);
%         REF = ampDataTA(k,ref_index);
%         REF = data(:,ref_index);

    
    window_size = 500;
    pixel_size = 80;
    ave_data = [];
    range = 1501:7001;
    for i = 1:1024
%        ave_data(i,:) = mean(data(i,specific_index-300:specific_index+300)); 
%        tmp = corrcoef(data(951,range), data(i,range));
%        ave_data(i) = tmp(1,2);
        ave_data(i,:) = (data(i, specific_index) - REF);
    end
    
    figure1=figure('Position', [100, 100, 100+window_size, 100+window_size]);
%     scatter(elecX(k),elecY(k),pixel_size,data(k,specific_index)-REF,'filled', 'o'); 
    scatter(elecX(k),elecY(k),pixel_size,ave_data(k),'filled', 'o'); 
    box on;
    ax=gca;ax.LineWidth=1.5;
    set(gca,'fontsize', 18);
        set(gca,'fontweight', 'bold');
%     S = std(data(k,specific_index)-REF);
	S = std(ave_data(k));

%     caxis([-1 1]);
    caxis([-3*S 3*S]); 
    h = colorbar;
%     ylabel(h, 'Slope (\muV    /ms)');
%     ylabel(h, 'f'''' (\muV/ms^2)');
    ylabel(h, 'Voltage (\muV)');
    set(gca,'FontSize',15)
    colormap('redblue');
%         colormap('jet');

    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
% ax_scl = 20;
%     axis([-ax_scl+1 ax_scl -ax_scl+1 ax_scl]);
    xlabel('mm');
    ylabel('mm');
    txt = ['Time = ' num2str((time)*1000, '%.3f ms')];
    title(txt);
    
    txt_filename = append('CS_RB_Nov13_2020_', num2str(t1), 's_',num2str((time)*1000, '%.0fms'), '.png');;
    print('-dpng','-r300',txt_filename);
    
    close(figure1);
    
    end
    %% Cross-Correlation - Plot at specific time
%         data = ampDataTA;
%         data = smoothdata(ampDataTA, 'gaussian', 50);
    
    for time_ms = [20:1:20]%[10:1:50, 60:10:100]
    ref_ms = 25.5;
    time = time_ms * 1e-3;
    specific_index = int32((time-timeStart)*fs);
    ref_index = int32((ref_ms*1e-3-timeStart)*fs);
%         REF = zeros(length(k),1);
%         REF = zeros(1024,1);
%         REF = ampDataTA(k,ref_index);
        REF = data(:,ref_index);

    
    window_size = 500;
    pixel_size = 80;
    ave_data = [];
    range = 1501:7001;
    for i = 1:1024
%        ave_data(i,:) = mean(data(i,specific_index-300:specific_index+300)); 
       tmp = corrcoef(data(951,range), data(i,range));
       ave_data(i) = tmp(1,2);
%         ave_data(i,:) = data(i, specific_index) - REF;
    end
    
    figure1=figure('Position', [100, 100, 100+window_size, 100+window_size]);
%     scatter(elecX(k),elecY(k),pixel_size,data(k,specific_index)-REF,'filled', 'o'); 
    scatter(elecX(k),elecY(k),pixel_size,ave_data(k),'filled', 'o'); 
    box on;
    ax=gca;ax.LineWidth=1.5;
    set(gca,'fontsize', 18);
        set(gca,'fontweight', 'bold');
%     S = std(data(k,specific_index)-REF);
	S = std(ave_data(k));

    caxis([0.5 1]);
%     caxis([-3*S 3*S]); 
    h = colorbar;
    ylabel(h, 'Correlation Coefficient');
    set(gca,'FontSize',15)
%     colormap('redblue');
        colormap('jet');

    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
% ax_scl = 20;
%     axis([-ax_scl+1 ax_scl -ax_scl+1 ax_scl]);
    xlabel('mm');
    ylabel('mm');
    txt = ['Time = ' num2str((time)*1000, '%.3f ms')];
%     title(txt);
    
    txt_filename = append('CCoeff2_Nov13_2020_', num2str(t1), 's_',num2str((time)*1000, '%.0fms'), '.png');;
    print('-dpng','-r300',txt_filename);
    
%     close(figure1);
    
    end
 %% Mapping plot of ampDataTaCMR
    fs = 2e4;
    raw_TA_data = downsample(ampDataTA', fs/1000)';
    fs = 1e3;
    range = (int32(-timeStart*fs)+5):(int32(-timeStart*fs)+50);
	x_size = 1000;
    y_size = 1000;
    line_width = 1;
    font_size = 10;
    figure1=figure('Position', [0, 0, x_size, y_size]);
    t_div = 50;
    v_div = 100;
    l_width = 0.6;
    
    raw_TA_rms = [];
    for i = 1:1024
        raw_TA_rms(i) = rms(raw_TA_data(i,range));
    end
    plot(raw_TA_rms(k))
    raw_TA_rms = raw_TA_rms/max(raw_TA_rms(k));
%     scatter(elecX(k),elecY(k),150,raw_TA_data(k,142),'filled', 's'); 
%     colorbar;
    
%     cm = colormap(jet); 
    c_index = (int32(-timeStart*fs)+28);
    max_amp = max(max(abs(raw_TA_data(k,c_index))));
    
    for i = k
        if raw_TA_data(i,c_index) > 0
        plot(elecX(i)-0.075+(1:length(range))/t_div, elecY(i)+raw_TA_data(i,range)/v_div, 'Color',[raw_TA_data(i,c_index)/max_amp,0,0], 'LineWidth', l_width);
        else
        plot(elecX(i)-0.075+(1:length(range))/t_div, elecY(i)+raw_TA_data(i,range)/v_div, 'Color',[0,0,-raw_TA_data(i,c_index)/max_amp], 'LineWidth', l_width);
        end
        hold on;
    end
    
    xlabel('Position (mm)');
    ylabel('Position (mm)');
    pbaspect([x_size y_size 1]);
    lim = 17   ;% - 0.075;
    axis([-lim lim -lim lim]);
    axis off;

    set(gca,'FontSize',font_size);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',line_width);
    
    txt = 'CS_mapping_Nov13_2020_N10';
    print('-dpng','-r600',txt);
    %% Plot at specific time for 1cm pitch, 3mm diameter
    time_ms = 43;
    time = time_ms * 1e-3;
    specific_index = int32((time-timeStart)*fs);
    
    num_c = 16;
    clinDataTA = zeros(num_c, 10001);
    clin_count = 1;
    num_points = 0;
    
    elecXclin = [];
    elecYclin = [];
    
        for j = -14.5:10:15.5
            for m = -14.5:10:15.5
                num_points = 0;
                for i = k
                    if ((elecX(i) <= j+1) && (elecX(i) >= j-1)) && ((elecY(i) <= m+1) && (elecY(i) >= m-1))
                        clinDataTA(clin_count,:) = clinDataTA(clin_count,:) + ampDataTA(i,:);
                        num_points = num_points + 1;
                    end
                end
                elecXclin(clin_count) = j;
                elecYclin(clin_count) = m;
                clinDataTA(clin_count,:) = clinDataTA(clin_count,:)/num_points;
                num_points
                clin_count = clin_count+1;
            end
        end
%         
%         for i = 1:16
%         plot(i*30+clinDataTA(i,:))
%         hold on;
%         end
%         elecXclin = elecXclin /1.5;
        
%         elecYclin = elecYclin /1.5;

    k_clin = 1:num_c;
    window_size = 500;
    pixel_size =1;
    figure1=figure('Position', [100, 100, 100+window_size, 100+window_size]);
    scatter(elecXclin(k_clin),elecYclin(k_clin),pixel_size,clinDataTA(k_clin,specific_index),'filled', 'o'); 
    box on;
    ax=gca;ax.LineWidth=1.5;
    set(gca,'fontsize', 18);
        set(gca,'fontweight', 'bold');


%     S = std(clinDataTA(k_clin,specific_index));
%     caxis([-5 5]);
    caxis([-3*S 3*S]); 
    h = colorbar;
    ylabel(h, 'Voltage (\muV)');
    set(gca,'FontSize',15)
    colormap('redblue');
    pbaspect([1 1 1]);
    ax_scl = 17;
    axis([-ax_scl ax_scl -ax_scl ax_scl]);
    ax = gca;               % get the current axis
ax.Clipping = 'on';    % turn clipping off
    xlabel('mm');
    ylabel('mm');
    txt = ['Time = ' num2str((time)*1000, '%.3f ms')];
    title(txt);
    %% Plot along y = const
%     vMargin = 1000;
    vOffset = 1.5;

    n_y = [];
    j = 1;
    y_pos = -3.5
    for i = 1:length(k)
        if (abs(elecY(k(i)) - y_pos)) < 0.1
            n_y(1, j) = elecX(k(i));
            n_y(2, j) = k(i);
            j = j + 1;
        end
    end
    j
%     Sorting n_y
    n_y = transpose(sortrows(transpose(n_y), 'descend'))
    
%     Plot along the line in n_y datas
    figure2 = figure('Position', [100, 0, 100+400, 0+800]);
    for i = 1:3:length(n_y(2,:))-2
    plot(((1:length(ampDataTA))/fs+timeStart)*1e3, vOffset*i+ampDataTA(n_y(2,i),:), 'linewidth', 2, 'color', [0 0 0]);
    hold on;
    end
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    xlim([-10 100]);
    ylim([-10 50]);

