% read_Intan_RHD2000_file_60Hz_notch
% AB_motor_v3_reader
% clear all
jNum = 1;
load('elecCoords_CS_32x32.mat');
elecX = elecCoords_CS_32x32(:,1);
elecY = elecCoords_CS_32x32(:,2);
% 329 - t_amplifier(1)
% load('k_vector_Linear_CS_11132020.mat');
load('k_vector_CS_11132020.mat');
% save('k_vector_CS_11132020.mat', 'k')
% load('CS_BPF.mat');
%% Referencing: Substracting minimum RMS 
% Common mode subtraction from minimum RMS channel

% 	for t1 = 1:2:59%52.59;%49.53;%4.7;%57.032%49.59+3%41%22;
%         t_amplifier(1)
%     t1 = 28.8;%79.3-60.032%49.6;
    motor_trial_num = 6;
                % 1     2     3     4     5     6      7    8    9
    motor_time = [298 303.4 309.6 313.5 320.9 325.76 335.1 350 313.16];
%     t1 = motor_time(motor_trial_num) - t_amplifier(1);
	t1 = 325.76 - t_amplifier(1);

    %motor_time(motor_trial_num)
    fs = 2e4;
    timeStart = -1;
    timeInterval = 2*abs(timeStart); % length(ampDataTA)/fs;%
    
    r1 = int32((t1+timeStart)*fs+1);
    r2 = int32((t1+abs(timeStart))*fs+1);
    range = r1:r2;

    ampData = amplifier_data(:,range);

%     ampData_RMS = [];
%     for i = 1:1024
%             ampData_RMS(i,1) = rms(ampData(i,:));
%             ampData_RMS(i,2) = i;
%             i
%     end
%     X = [ampData_RMS(k,1) ampData_RMS(k,2)];
%     Y = sortrows(X, 'ascend');
    
    num_CM = length(k); % Number of minimum RMS channels to take common mode
    ampDataCM = zeros(1024,length(range));
    
    for i = 1:num_CM
    ampDataCM = ampDataCM + ampData(k(i),:);%ampData(Y(i, 2),:);
    display(sprintf('Common mode summation in progress.. %.1f%% completed.', 100*i/num_CM));
    end
    ampDataCM = ampDataCM/num_CM;

    ampDataCMR = [];
    ampDataCMR = ampData - ampDataCM;

    %% HG filtering
    fs = 2e4;
    N = 4;
    lf = [1 4 8 13 30 70 55 100 9 10]% 190 250 500 750]; % lower bound in Hz
    hf = [4 8 12 30 70 190 95 190 18 45]% 250 500 750 2000]; % Upper bound in Hz
    
    i_freq = 4;% 5: Low gamma / 6: High Gamma Activity
    d_B = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
    'SampleRate',fs);

    i_freq = 6;% High Gamma Activity
    d_HG = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
    'SampleRate',fs);
    
    % High gamma filtering
    data_ds = [];
    data_HG_ds = [];
    data_HG_amp_ds = [];
    data_B_amp_ds = []
    data_B_ds = [];
    for i = 1:1024
    data_HG_temp = filtfilt(d_HG, ampDataCMR(i,:));
    data_B_temp = filtfilt(d_B, ampDataCMR(i,:));
	data_HG_ds(i,:) = downsample(data_HG_temp, fs/1e3);
    data_ds(i,:) = downsample(ampDataCMR(i,:), fs/1e3);
    data_HG_amp_ds(i,:) = abs(hilbert(downsample(data_HG_temp, fs/1e3)));
    data_B_amp_ds(i,:) = abs(hilbert(downsample(data_B_temp, fs/1e3)));
    data_B_ds(i,:) = downsample(data_B_temp, fs/1e3);
    display(sprintf('Hilbert transform in progress.. %.1f%% completed.', 100*i/1024));
    end
    
%       save(append('Motor_11Nov2020_Grab_trial5_at325p76', '.mat'), 'data_ds', 'data_HG_amp_ds', 'data_B_amp_ds', 'data_B_ds');
%       save(append('Motor_11Nov2020_Background', '.mat'), 'data_ds_BG', 'data_HG_amp_ds_BG', 'data_B_amp_ds_BG', 'data_B_ds_BG');

%     for i = 1:100:1024
%     plot((1:length(data_HG_ds))/fs,0.2*i+data_HG_ds(i,:))
%     hold on;
%     end
 %% Trial average
% data_ds_TA = zeros(1024, 2001);
% data_HG_amp_ds_TA = zeros(1024, 2001);
% data_B_amp_ds_TA = zeros(1024, 2001);
% data_B_ds_TA = zeros(1024, 2001);
% 
%  for motor_trial_num = 1:6
%  load(append('Motor_11Nov2020_Grab_trial', num2str(motor_trial_num), '.mat'));
%  data_ds_TA = data_ds_TA + data_ds;
%  data_HG_amp_ds_TA = data_HG_amp_ds_TA + data_HG_amp_ds;
%  data_B_amp_ds_TA = data_B_amp_ds_TA + data_B_amp_ds;
%  data_B_ds_TA = data_B_ds_TA + data_B_ds;
% 
%  motor_trial_num
%  end
%  data_ds_TA = data_ds_TA/6;
%  data_HG_amp_ds_TA = data_HG_amp_ds_TA/6;
%  data_B_amp_ds_TA = data_B_amp_ds_TA/6;
%  data_B_ds_TA = data_B_ds_TA/6;
 
 %Substracting background
 load('Motor_11Nov2020_background.mat');
%  data_ds_TA = data_ds_TA - data_ds;
%  data_HG_amp_ds_TA = data_HG_amp_ds_TA - data_HG_amp_ds;
%  data_B_amp_ds_TA = data_B_amp_ds_TA - data_B_amp_ds;
%  data_B_ds_TA = data_B_ds_TA - data_B_ds;
%% Filtering the 1kHz DownSampled Data
%     fs = 1e3;
%     N = 4;
%     lf = [1 4 8 13 30 70 55 100 9]% 190 250 500 750]; % lower bound in Hz
%     hf = [4 8 12 30 70 190 95 190 18]% 250 500 750 2000]; % Upper bound in Hz
%     
%     i_freq = 4;% 5: Low gamma / 6: High Gamma Activity
%     d_B = designfilt('bandpassiir','FilterOrder',N, ...
%     'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
%     'SampleRate',fs);
% 
%     i_freq = 5;%  5: Low gamma 
%     d_LG = designfilt('bandpassiir','FilterOrder',N, ...
%     'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
%     'SampleRate',fs);
% 
%     i_freq = 6;% High Gamma Activity
%     d_HG = designfilt('bandpassiir','FilterOrder',N, ...
%     'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
%     'SampleRate',fs);
%     data_LG_amp_ds = [];
%     
%     for i = 1:1024
%     data_HG_temp = filtfilt(d_HG, data_ds(i,:));
%     data_B_temp = filtfilt(d_B, data_ds(i,:));
%         data_LG_temp = filtfilt(d_LG, data_ds(i,:));
% 
% % 	data_HG(i,:) = abs(hilbert(data_HG_temp));
% %     data_ds(i,:) = downsample(ampDataCMR(i,:), fs/1e3);
%     data_HG_ds(i,:) = data_HG_temp;%, 'gaussian', 20);
%     data_HG_amp_ds(i,:) = abs(hilbert(data_HG_temp));%, 'gaussian', 20);
%     data_B_amp_ds(i,:) = abs(hilbert(data_B_temp));
%     data_B_ds(i,:) = data_B_temp;
%     data_LG_amp_ds(i,:) = data_LG_temp;
%     display(sprintf('Hilbert transform in progress.. %.1f%% completed.', 100*i/1024));
%     end
    %% Sort channels and Plot
    offset = 5;
    fs = 1e3;
    data = data_HG_amp_ds;
    
    pos_xy = 10;
    window_size_x = 1900;
    window_size_y = 1000;
      figure1=figure('Position', [pos_xy, pos_xy, pos_xy+window_size_x, pos_xy+window_size_y]);
%   k_sorted(19,2)= 640
%	k_sorted(18,2) = 630 Motor
%     k_sorted(3,2) = 476 Sensory
%     count = 0;
%     for i_num = [630, 476]
%     plot((1:length(data))/fs+timeStart, count*10 + data(i_num,:) )
%     count = count + 1;
%     hold on;
%     end
    
    for x_pos = [15.5]%-15.5:1:15.5;
        
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
    
    figure1=figure('Position', [pos_xy, pos_xy, pos_xy+window_size_x, pos_xy+window_size_y]);
    
    
    for i = 1:length(k_temp)
    plot((1:length(data))/fs+timeStart,offset*i+ data(k_sorted(i,2),:) )
    hold on;
    end
    
    i=1
	plot((1:numPoints)/408.16 + 221.4 + 35.7900-(t_amplifier(1)+t1), 5e2*transpose(adcVals(i,:)-adcVals(i,1))-50);
    hold on;
    i=2:6
	plot((1:numPoints)/408.16 + 221.4 + 35.7900-(t_amplifier(1)+t1), 1e2*transpose(adcVals(i,:)-adcVals(i,1))-50);
    xlim([-0.5 0.5]);
    
    xlim([timeStart -timeStart]);
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ylim([-offset-100 offset*(length(k_temp)+1)]);
    txt = ['X = ' num2str(x_pos, '%.1f mm') '_' num2str(t_amplifier(1)+t1), 's_'];
    title(txt);
    
    txt_filename = append('Sensory_shift_Beta_Nov13_2020_', num2str(t1), 's_', num2str(x_pos), 'mm_', '.png');;
%     print('-dpng','-r600',txt_filename);

    pause(2) %in seconds
    hold on;

%         close(figure1);

    end
    
    %% Plot the Tactile Neuron
%                 i = 1;
        data = data_HG_amp_ds;
        fs = 1e3;
            plot(t_amplifier(1)+t1+timeStart+(1:length(data(923,:)))/fs, data(923,:))
            plot(timeStart+(1:length(data(923,:)))/fs, data(923,:))
                        plot(data(923,:))

            xlabel('Time (s)');
            ylabel('Voltage (\muV)');
%             xlim([-10 10]);
    %% Plot RMS map
        M = 10;
    f_size = 15;
    t_center = -1%-0.1; % s
	ref_time = -1%0.2;
    RMS_calc_time = 0.1;
	Ref_RMS_calc_time = 0.9;%0.1;
	dtime = 0.05;
    fs = 1e3;
    c_max = 3;
    c_min = 1;


	filename = append('Motor_HG_Trial2_325p76_BG1s3_RB_Nov13_2020_', num2str(t_amplifier(1)+t1+t_center), 's_', num2str(jNum), '.gif');;
    jNum = jNum + 1;
        data = data_HG_ds;
%         data_ref = data_HG_amp_ds_BG;

    window_size_x = 800;
	window_size_y = window_size_x;
    pixel_size = 150;%40*power(window_size_x/400,2);
    figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);
        set(gcf,'color','w');

    data_sigma = zeros(1024, 2001);
    
    for i_gif = 1:1:(1900)
%     i_gif
    plus = i_gif*0.001; %ms
    timeSt = t_center + plus;
    
	range = int32((-timeStart+timeSt)*fs):int32((-timeStart+timeSt+RMS_calc_time)*fs);
    range_ref = int32((-timeStart+ref_time)*fs)+1:int32((-timeStart+ref_time+Ref_RMS_calc_time)*fs);
%     range_ref(end)
    data_rms = [];
    for i = 1:1024
                data_rms(i,1) = rms(data(i, range))/rms(data(i, range_ref));
                data_sigma(i, i_gif) = data_rms(i,1);
%         data_rms(i,1) = rms(data(i, range))/rms(data(i, range_ref));
%         data_rms(i,1) = rms(abs(hilbert(data(i, range))))/rms(abs(hilbert(data(i, range_ref))));
%        data_rms(i,1) = rms(abs(hilbert(data(i, range)))) - rms(abs(hilbert(data(i, range_ref))));% - rms(data_ref(i, range_ref));
    end
    
    scatter(elecX(k),elecY(k),pixel_size,data_rms(k),'filled', 'o'); 
%     grid on;
    ax = gca;
    ax.GridColor = [0 0 0];
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.8;
    
    box on;
    ax=gca;ax.LineWidth=1;
    set(gca,'fontsize', f_size);
    set(gca,'fontweight', 'bold');
    S = std(data_rms(k));
%             caxis([0 5*S]);
%     M = 2;%max(data_rms(k));
    caxis([c_min c_max]);
%     caxis([0 M]);
%     caxis([-3*S 3*S]); 
    h = colorbar;
	ylabel(h, 'sigma (\sigma)');
%     ylabel(h, 'Voltage (\muV)');
    set(h,'FontSize',f_size)
    set(gca,'FontSize',f_size)
%     colormap('jet');
%     colormap(flipud(hot))
%     colormap('redblue');
%     ZtoO = power((0:127)/127, 3);
%     OtoZ = ones(1,128) - ZtoO;
	R = power((0:255)/255, 3);
    myColorMap = [ones(256,1), ones(256,1)-R', ones(256,1)-R'];
% %     myColorMap = [[OtoZ linspace(0,0,128)]', [OtoZ linspace(0,0,128)]', [linspace(1,1,128) OtoZ]'];
    colormap(gca, myColorMap);
    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
    xlabel('mm');
    ylabel('mm');
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'Color','k')

%     caxis([0 10]);
    txt = [num2str((plus+t_center)*1e3, '%.0f ms')];
%     txt = [num2str((plus+t_center)*1e3, '%.0f ms')];
%     txt = [num2str(t_amplifier(1)+t1+timeSt, '%.2f s / ') num2str((plus+t_center)*1e3, '%.0f ms')  ' / 0~' num2str((M), '%.0f') '\muV'];
    title(txt);        
    
%     if rem(i_gif, 1) == 0
%          print('-dtiff','-r300',append('Motor_HG_RB_map_BG1s3_trial5_', num2str((plus+t_center)*1e3, '%.0f ms')));
%     end


%     txt_filename = append('Motor_HGA_Nov13_2020_', num2str(t_amplifier(1)+t1+timeSt), 's', '.png');;

     % Capture the plot as an image 
      frame = getframe(figure1); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
      if i_gif == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime); 
      end 
        hold off;
    
    end
    %     print('-dpng','-r300',txt_filename);
    close(figure1);
    %% Motion plot
            window_size_x = 1900;
    window_size_y = 300;
	figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);

    
            i=3
	plot((1:numPoints)/408.16 + 221.4 + 35.7900-(t_amplifier(1)+t1), smoothdata(1e2*transpose(adcVals(i,:)-adcVals(i,1)), 'gaussian', 10),...
        'linewidth', 5, 'color', 'k');
    xlim([-0.5 0.5]);
            set(gca,'YTick',[])
    set(gca,'XTick',[])
        print('-dpng','-r300','Motion_1s.png');

%% ImageSC
%     AB_motor_v3_reader_Motor_Nov13_2020
    app_txt = 'Trial2_375p76s'
	window_size_x = 2400;
    window_size_y = 300;
	figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);

    imagesc(data_imagesc)
	colormap(redblue);
    colorbar;
    caxis([0 10]);
    
%     time1 = -0.1-0.5;
    time1 = -0.5;
    time2 = time1+1;
	xlim([(-timeStart+time1)*Fs (-timeStart+time2)*Fs]);
	axis off;
%     print('-dpng','-r300',append('Beta_amp_', app_txt,'_13to30Hz_1s_Nov13_2020.png'));
        
	window_size_x = 2400;
    window_size_y = 300;
	figure2=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);

    imagesc(data_imagesc_B)
	colormap(redblue);
    colorbar;
    caxis([-50 50]);
	xlim([(-timeStart+time1)*Fs (-timeStart+time2)*Fs]);
% 	axis off;
%     print('-dpng','-r300',append('HG_amp_', app_txt,'_70to190Hz_1s_Nov13_2020.png'));

%     xlim([5e2 15e2]);
%     print('-dpng','-r300','HG_amp_map_300to326s_70to190Hz_1s_Nov13_2020.png');

        %% Plot Along a Line
    i_num = 640;
    Fs = 1e3;
        window_size_x = 2400;
    window_size_y = 1200;
	figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);

%     plot(data_B_amp_ds(i_num, :));
%     hold on;
    data = data_HG_amp_ds;
    data_B = data_B_ds;
%     data = abs(hilbert(data_B_ds));
    data_imagesc = [];
	data_imagesc_B = [];


    v_offset =  -30;
        count = 0;

    for i_scan = 1:32
%         x_c = 15.5 - i_scan;
        x_c = 1.5;% - i_scan;

        y_c = 15.5 - i_scan;
    R = 0.1;
    

    for i = 1:length(k)
        if ((elecX(k(i)) > x_c-R) && (elecX(k(i)) < x_c+R)) ...
                && ((elecY(k(i)) > y_c-R) && (elecY(k(i)) < y_c+R))
%             k(i)

            data_imagesc = [data_imagesc data(k(i),:)'];
            data_imagesc_B = [data_imagesc_B data_B(k(i),:)'];

            plot(((1:timeInterval*Fs+1)+timeStart*Fs)/Fs, count*v_offset + data_B(k(i),:), 'color', 'k', 'linewidth', 1);
            count = count + 1;
            hold on;
            
        end
    end
        i_scan
    
    end
    data_imagesc = data_imagesc';
    data_imagesc_B = data_imagesc_B';


    
%         xlim([-5e2 5e2]);

    t_plot_center = 0;%85;%420;%-255;
    t_off = 100;
%     xlim([t_plot_center-t_off t_plot_center+t_off]);
    title(append(num2str(t_plot_center, '%.2f'), ' ms'));
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
%     ylim([-400 50]);
%     set(gca,'YTick',[])
%     txt_filename = append('Beta_waveforms2_', num2str(t_plot_center, '%.0f'),'ms_', num2str(t_amplifier(1)+t1+timeSt), 's', '.png');;
%     print('-dpng','-r300',txt_filename);

    
    %% Plot averaged waveforms
    i_num = 640;
%     plot(data_B_amp_ds(i_num, :));
%     hold on;

    x_c = 12.5;
    y_c = 3.5;
    R = 0.1;
    data_HG_ave = zeros(1, length(data_HG_amp_ds(i,:)));
        data_B_ave = zeros(1, length(data_B_amp_ds(i,:)));
        data_ave = zeros(1, length(data_ds(i,:)));

    count = 0;
    for i = 1:length(k)
        if ((elecX(k(i)) > x_c-R) && (elecX(k(i)) < x_c+R)) ...
                && ((elecY(k(i)) > y_c-R) && (elecY(k(i)) < y_c+R))
%             k(i)
            data_HG_ave = data_HG_ave + data_HG_amp_ds(k(i),:);
            data_B_ave = data_B_ave + data_B_amp_ds(k(i),:);
            data_ave = data_ave + data_ds(k(i),:);

            count = count + 1;
        end
    end
    data_HG_ave = data_HG_ave / count;
	data_B_ave = data_B_ave / count;


	plot(data_HG_ave(:)*10);
    hold on;
	plot(data_B_ave(:));
        hold on;
	plot(data_ave(:));

%%     spectrogram(data_B_ave)
fs = 1e3;
spectrogram(data_ave(800:1200),11,10,5e2,fs, 'yaxis');
caxis([-15 35]);
    
 %% Plot R2 map
    window_size_x = 900;
	window_size_y = 900;
	figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);

    timeSt = -6.85
    ref_time = -1;
    RMS_calc_time = 0.1;
	Ref_RMS_calc_time = 0.1;
    
    data = data_HG_r_ds;
    fs = 1e3;

	range = int32((-timeStart+timeSt)*fs):int32((-timeStart+timeSt+RMS_calc_time)*fs);
    range_ref = int32((-timeStart+timeSt+ref_time)*fs):int32((-timeStart+timeSt+Ref_RMS_calc_time+ref_time)*fs);
    length(data(i,range))
    length(data(i,range_ref))
    
%     i = 923
%     hold on;
%     plot(data(i, range_ref))
%     
    pixel_size = 150;
   
    
    data_R2 = [];
    for i = 1:1024
%        corr_temp = cov(data(i, range_ref), data(i, range));
%        std_prod = std(data(i, range_ref))*std(data(i, range));
%        data_R2(i,1) = power(corr_temp(1,2)/std_prod, 2);
       corr_temp = corrcoef(data(i, range_ref), data(i, range));
       data_R2(i,1) = -log10(power(corr_temp(1,2), 2));
       
%        d1 = [data(i, range_ref) data(i, range)];
%     l1 = length(data(i, range_ref));
%     ave1 = mean(data(i, range_ref));
%     ave2 = mean(data(i, range));
%     d2 = [ ave1*ones(1,l1) ave2*ones(1,l1)];
%     corr_temp = corrcoef(d1, d2);
% 	data_R2(i,1) = power(corr_temp(1,2), 2);
     
%         data_R2(i,1) = var(data(i, range)) - var(data(i, range_ref));
       display(sprintf('Correlation Coeff. calc. in progress.. %.1f%% completed.', 100*i/1024));

    end
    
    scatter(elecX(k),elecY(k),pixel_size,data_R2(k),'filled', 'o'); 
    grid on;
    ax = gca;
    ax.GridColor = [0 0 0];
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.8;
    
    box on;
    ax=gca;ax.LineWidth=1;
    set(gca,'fontsize', 15);
    set(gca,'fontweight', 'bold');
%     S = std(data_R2(k));
%             caxis([0 5*S]);
    M = max(data_R2(k));
%     caxis([0 M]);
%     caxis([-3*S 3*S]); 
    h = colorbar;
    ylabel(h, 'R^2');
    set(gca,'FontSize',15)
    colormap(flipud(hot))
%     colormap('redblue');
    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
    xlabel('mm');
    ylabel('mm');
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'Color','k')

%     caxis([0 10]);

    txt = [num2str(t_amplifier(1)+t1+timeSt, '%.2f s') ' / 0~' num2str((M), '%.0f') '\muV'];
    title(txt);        

    txt_filename = append('Sensory_HG_R2_Nov13_2020_', num2str(t_amplifier(1)+t1+timeSt), 's', '.png');;
    print('-dpng','-r300',txt_filename);
%     close(figure1);
   
    %% Plot at specific time
    
    timeSt = -7.19
    
    for it_loop = 0:10
        
    data = data_HG_amp_ds;
    fs = 1e3;
    
    center_time = (timeSt + 0.02*it_loop )*1e3;%pkIndex(i_pk)-5;
    window_size_x = 1900;
	window_size_y = 900;

    pixel_size = 13;
    figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);
    
    for i_loop = 1:18%, 60:10:100]1528 -595 1073
%     ref_ms = 5;
    i_ms = i_loop * 1;
    time_ms = center_time + i_ms;
    time = time_ms * 1e-3;
    specific_index = int32((time-timeStart)*fs);
%     ref_index = int32((ref_ms*1e-3-timeStart)*fs);
%         REF = zeros(length(k),1);
%         REF = ampDataTA(k,ref_index);
    
    
    subplot(3,6,i_loop);
    scatter(elecX(k),elecY(k),pixel_size,data(k,specific_index),'filled', 'o'); 
    
    box on;
    ax=gca;ax.LineWidth=1;
    set(gca,'fontsize', 10);
    set(gca,'fontweight', 'bold');
    S = std(data(k,specific_index));
%     caxis([0 15]);
%     caxis([-3*S 3*S]); 
%     h = colorbar;
%     ylabel(h, 'Voltage (\muV)');
%     set(gca,'FontSize',10)
    colormap jet;
%     colormap('redblue');
    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
%     xlabel('mm');
%     ylabel('mm');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
        caxis([0 3*S]);

    txt = ['+' num2str(i_ms, '%.0f ms') ' / 0~' num2str((S), '%.0f') 'mV'];
    title(txt);    
    end    
    
    txt = ['Time = ' num2str((t_amplifier(1)+t1+center_time*1e-3), '%.3f')  's'];
    sgtitle(txt)
    

    txt_filename = append('Sensory_HG_Nov13_2020_', num2str(t_amplifier(1)+t1+center_time*1e-3), 's', '.png');;
    print('-dpng','-r300',txt_filename);
    close(figure1);

    end
    
%     end

%% Spectogram
i_from_bottom = 20
% 
% plot(data_B_amp_ds(k_sorted(i_from_bottom,2),range))

figure5=figure('Position', [50, 50, 50+1900, 50+500]);
data = data_ds;%ampDataCMR;%
fs = 1e3;
t_spec = -0.2;
t_intv = 0.5;
range = int32((-timeStart+t_spec)*fs+1):int32((-timeStart+t_spec+t_intv)*fs+1);
% plot(data(630,range)) 476
spectrogram(data(476,range),21,20,5e2,fs, 'yaxis');
caxis([-10 30]);
% ylim([1e-3 1])
ylim([1 200])
%     colormap(flipud(hot))

colormap jet

%     plot((1:length(data))/fs+timeStart,offset*i+ data(k_sorted(i,2),:) )
ax = gca;
% ax.YScale = 'log';

% yt = get(gca, 'YTick');
% ax.YTick = yt;
% ax.YTickLabel = yt*1E3;
xt = get(gca, 'XTick');
ax.XTick = xt;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% A = -200*ones(1, length(xt));
% ax.XTickLabel = xt + A;
h = colorbar;
ylabel(h, 'dB/Hz');
set(gca,'FontSize',15);
    txt = append('Spectrogram_InvHot_shift_', filename, '_t_', num2str(t_amplifier(1)+t1-timeStart+t_spec, '%.0f'), '_ch_', num2str(k_sorted(i_from_bottom,2), '%.0f'), '.png');
%     saveas(gcf,txt);
%     close(figure5);
% 
