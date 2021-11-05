
   jNum = 100;
    
     %% Load Method 2 data
    fs = 2e4;
%     load('Depth_Stim_Method2_May14_2021');
    t_ini = 0.0592-(50)*1e-3;%2.001-50*1e-3;%3.845-50*1e-3;
    range = int32(t_ini*fs):int32((t_ini+0.5)*fs)
    d = ampDataCMR';

    % Load Mapping info & k-vector
    ch_num = 768;
    load('elecCoords_v14_32x32');
    elecX = elecCoords_v14_32x32(:,1);
    elecY = elecCoords_v14_32x32(:,2);
    load('k_vector_surf_v14_May14_2021.mat');
    
   
    %% Band pass filtering
    fs = 2e4;
    Fs = 1e3;
    N = 4;
    lf = [1 4 8 13 30 70 55 100 9 10]% 190 250 500 750]; % lower bound in Hz
    hf = [4 8 12 30 70 190 95 190 18 59]% 250 500 750 2000]; % Upper bound in Hz
    
    i_freq = 10;% 5: Low gamma / 6: High Gamma Activity
    d_B = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf(i_freq),'HalfPowerFrequency2',hf(i_freq), ...
    'SampleRate',fs);

    data_B_ds = [];
	data_B_amp_ds = [];

    for i = 1:ch_num
    data_B_temp = filtfilt(d_B, d(range, i)');
	data_B_amp_temp = abs(hilbert(data_B_temp));
    data_B_ds(i,:) = downsample(data_B_temp, fs/Fs);
    data_B_amp_ds(i,:) = downsample(data_B_amp_temp, fs/Fs);
    display(sprintf('Hilbert transform in progress.. %.1f%% completed.', 100*i/ch_num));
    end

    %% Convert to 32 x 32
    % Feeding the data
    data = data_B_ds;%downsample(xB', 20000/1000)';;%downsample(xB', 20000/1000)';%artiData;%data_Mc_HG_Processed;data_Mc_HG_Processed;%
    data_amp = data_B_amp_ds;
    % Data length
    d_length = int32(size(data));
    d_length = d_length(2);
    
    % Interpolation
    interp_x_segments = 32; % Number of segments in x-axis
    interp_y_segments = interp_x_segments ; % Number of segments in y-axis
    % set up new X, Y, and Z coordinates for interpolated output
    [xmin,xmax] = bounds(elecX);
    [ymin,ymax] = bounds(elecY);
    [xGrid,yGrid] = meshgrid(linspace(xmin,xmax,interp_x_segments),linspace(ymin,ymax,interp_y_segments));

    zGrid = zeros(32,32, d_length);
	zGrid_amp = zeros(32,32, d_length);

    for t_loop = 1:d_length
    Finter= scatteredInterpolant(elecX(k),elecY(k),data(k, t_loop));
    zGrid(:,:,t_loop) = Finter(xGrid,yGrid);
    
    Finter= scatteredInterpolant(elecX(k),elecY(k),data_amp(k, t_loop));
    zGrid_amp(:,:,t_loop) = Finter(xGrid,yGrid);
    t_loop
    end
    %% Phase Gradient Vector field plot
    addpath('wave')
    addpath('wave\analysis')
    addpath('wave\plotting')

	time_ini = 1; % in ms

    % parameters
    Fs = 1000;%Hz
    image_size = 16; %px
    pixel_spacing = 1.5; %a.u.
    direction = +1; % +1/-1

    % generate data
    xf = zGrid;
    % z-score data
    xf = zscore_independent( xf );
    % form analytic signal
    xph = analytic_signal( xf );
    % calculate instantaneous frequency 
    [wt,signIF] = instantaneous_frequency( xph, Fs );
    % calculate phase gradient
    [pm,pd,dx,dy] = phase_gradient_complex_multiplication( xph, pixel_spacing, signIF );
  
    pd_unwrap = [];
    for i=1:32
        for j=1:32
    pd_unwrap(i,j,1:d_length) = reshape(unwrap(reshape(pd(i,j,1:d_length), [1 d_length])), [1 1 d_length]);
    pd_unwrap(i,j,1:d_length) = pd_unwrap(i,j,1:d_length) - pd_unwrap(i,j,time_ini) + pd(i,j,time_ini);
        end
        i*j
    end
    
    %% RMS video GIF plotting
    jNum = jNum + 1; 
%     custFolderName = 'C:\Beta_wave_Sep25_2020\';
    fName = 'Baseline_Propagating_wave_1to60Hz_14May2021';
    v_max = 300;
    v_min = -v_max;
    
    t_window = [1:2:301]
    figure2 = figure('Position', [0, 0, 0+800, 0+800]);
    set(gcf,'color','w');

    filename = append(fName, '_' ,num2str(jNum), '.gif');
for  i = 1:length(t_window)
    time_stamp = t_window(i);
    [C,h]= contourf(1:32,1:32,zGrid(:,:,time_stamp),100);
    set(h,'LineColor','none');
%     colormap('redblue');
%     colormap(flipud(hot));
%     R = power((0:255)/255, 1)
%     myColorMap = [ones(256,1), ones(256,1), ones(256,1)-R'];
%     colormap(gca, myColorMap);
    colormap(gca, redblue);
    axis([1 32 1 32]);
    caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    
    txtTitle = append(fName, '_at_', num2str(time_stamp
    ), '_ms' );
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([32 32 1]);
    set(gca,'YTick',0:4:31);
    set(gca,'YTickLabel',((0:4:31)-16)*1.5);
    set(gca,'XTick',0:4:31);
    set(gca,'XTickLabel',((0:4:31)-16)*1.5);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',1);
    set(h,'linewidth',1);
    set(h,'FontSize',10);
    set(h,'FontWeight','bold');
    set(ht,'FontWeight','bold');

    hold on;
    % plot resulting vector field
    plot_vector_field_14May2021( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
    drawnow 
    
%         % Save as image
%         h = gcf;
%         set(h,'PaperOrientation','landscape');
        txt = append(txtTitle, '.png');
%     if rem(i, 1) == 0
%          print('-dpng','-r300',txtTitle)
%     end

      % Capture the plot as an image 
      frame = getframe(figure2); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      dtime = 0.15;
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime); 
      end 
        hold off;

        
        end
        close(figure2)
    %% Scatter Plot  video GIF plotting
    jNum = jNum + 1; 
%     custFolderName = 'C:\Beta_wave_Sep25_2020\';
    fName = 'Amp_mapping_25Sep2020_Baseline_10to60Hz';
    v_max = 100;
    v_min = 0;
    lim = 17;
	fs = 1e3;
    xy_offset = 500;
    pixel_size = 15;
    
	dtime = 0.05;
    f_size = 8;


    
    t_window = [1:10:4001]
    figure2 = figure('Position', [xy_offset, xy_offset, 300, 300]);
    set(gcf,'color','w');

    filename = append(fName, '_' ,num2str(jNum), '.gif');
for  i = 1:length(t_window)
    time_stamp = t_window(i);
    
    scatter(elecX(k),elecY(k), pixel_size, data_amp(k, time_stamp),'filled', 'o'); 
%     colormap('redblue');
    colormap(flipud(hot));
%     R = power((0:255)/255, 1)
%     myColorMap = [ones(256,1), ones(256,1), ones(256,1)-R'];
%     colormap(gca, myColorMap);
    axis([-lim lim -lim lim]);
    caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    
	txtTitle = append(num2str((time_stamp)/fs, '%.2f'), ' s' );
%     txtTitle = append(fName, '_', num2str(jNum), '_at_', num2str(time_stamp), '_ms' );
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([32 32 1]);

    set(gca,'FontSize',f_size);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',1);
    set(h,'linewidth',1);
    set(h,'FontSize',f_size);
    set(h,'FontWeight','bold');
    set(ht,'FontWeight','bold');
    box on;

    hold on;
    
%         % Save as image
%         h = gcf;
%         set(h,'PaperOrientation','landscape');
%         txt = append(txtTitle, '.png');
%     if rem(i, 1) == 0
%          print('-dpng','-r300',txtTitle)
%     end

      % Capture the plot as an image 
      frame = getframe(figure2); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime); 
      end 
        hold off;

        
        end
        close(figure2)
    