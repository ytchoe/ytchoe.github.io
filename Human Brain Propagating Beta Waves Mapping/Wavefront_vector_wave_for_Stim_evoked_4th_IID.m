    % Baseline - Spontaneous activities
%   save('spontaneous_epileptiform_data.mat', 'data');
%     load('spontaneous_epileptiform_data.mat');
    load('data2_stim_evoked_4th_epilepsy_Sep25_2020.mat');
    data = data_scatter_amp';

    % Electrode mapping information (X, Y)    
    load('elecCoords_CS_32x32.mat');
    elecX = elecCoords_CS_32x32(:,1);
    elecY = elecCoords_CS_32x32(:,2);
    % k vector contains the channel numbers with impedance below 100kOhm @ 1kHz
    load('k_vector_Sep25_2020.mat');
    %% Band pass filtering
    fs = 1e3; % Sampling rate
    N = 4; % Filter order
    lf = 10; % lower bound in Hz
    hf = 59; % Upper bound in Hz
    
    % Filter: Butterworth 4th order filter
    d_B = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',lf,'HalfPowerFrequency2',hf, ...
    'SampleRate',fs);

    data_B_ds = []; % Band-pass filtered data
	data_B_amp_ds = []; % Amplitude

    for i = 1:1024
    data_B_ds(i,:) = filtfilt(d_B, data(:, i)'); % Band-pass filtered data
	data_B_amp_ds(i,:) = abs(hilbert(data_B_ds(i,:))); % Amplitude
    
    display(sprintf('Filtering in progress.. %.1f%% completed.', 100*i/1024));
    end

    %% Interpolation, Convert to 32 x 32 interpolated array
    % Data length
    d_length = int32(length(data_B_ds));
    
    num_interp = 32;
    % Interpolation
    interp_x_segments = num_interp; % Number of segments in x-axis
    interp_y_segments = interp_x_segments ; % Number of segments in y-axis
    % set up new X, Y, and Z coordinates for interpolated output
    [xmin,xmax] = bounds(elecX);
    [ymin,ymax] = bounds(elecY);
    [xGrid,yGrid] = meshgrid(linspace(xmin,xmax,interp_x_segments),linspace(ymin,ymax,interp_y_segments));

    zGrid = zeros(num_interp,num_interp, d_length);
	zGrid_amp = zeros(num_interp,num_interp, d_length);
    zGrid_amp_smooth = zeros(num_interp,num_interp, d_length);
    
    for t_loop = 1:d_length
    Finter = scatteredInterpolant(elecX(k),elecY(k),data_B_ds(k, t_loop));
    zGrid(:,:,t_loop) = Finter(xGrid,yGrid); % interpolated array
    display(sprintf('Interpolation in progress.. %.1f%% completed.', 100*t_loop/d_length));
    end
    %% Phase Gradient Vector field plot
    addpath('wave')
    addpath('wave\analysis')
    addpath('wave\plotting')

    % parameters
    fs = 1000; % Sampling rate in Hz
    time_ini = 1; % in ms
    image_size = 16; % px
    pixel_spacing = 1; % a.u.
    direction = +1; % +1/-1

    % generate data
    xf = zGrid;
    % z-score data
    xf = zscore_independent( xf );
    % form analytic signal
    xph = analytic_signal( xf );
    % calculate instantaneous frequency 
    [wt,signIF] = instantaneous_frequency( xph, fs );
    % calculate phase gradient
    [pm,pd,dx,dy] = phase_gradient_complex_multiplication( xph, pixel_spacing, signIF );
    display(sprintf('Phase gradient calculation completed.'));

    %% RMS video GIF plotting
    txtTitle2 = 'Stim_Evoked_IID_4th_streamlines2';
    fName = 'Wave_front'; % File name
    v_max = 100; % Max voltage plot range
    v_min = 0; % min voltage plot range
    pixel_size = 150; % Scatter plot pixel size
    frame_rate = 10; % GIF frame rate
    
    t_interval = 100; % Time interval between the frame in 'ms'
    t_window = [1:1:8001]; % Time window in 'ms'
    time_zero = 2500;
    figure2 = figure('Position', [0, 0, 0+800, 0+800]);
    set(gcf,'color','w');

%     filename = append(fName, '.gif');

for  i = 2675%325%00:10:2400%1:length(t_window)
    time_stamp = time_zero + t_window(i); % time in 'ms'
    % Scatter plot
% 	scatter(elecX(k)+16.5,elecY(k)+16.5, pixel_size, data_B_amp_ds(k, time_stamp),'filled', 'o'); 
    colormap(flipud(hot));
    
    % Scatter plot: Settings
    caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    txtTitle = append(fName, '_at_', num2str(time_stamp-time_zero), '_ms' );
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([1 1 1]);
    set(gca,'YTick',0:4:32);
    set(gca,'YTickLabel',((0:4:32)-16));
    set(gca,'XTick',0:4:32);
    set(gca,'XTickLabel',((0:4:32)-16));
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',1);
% 	axis([-1 34 -1 34]);
    box on;
    % Colorbar: setting
    set(h,'linewidth',1);
    set(h,'FontSize',10);
    set(h,'FontWeight','bold');
    set(ht,'FontWeight','bold');

    hold on;

    % Plot Vector Field
    ph = exp( 1i .* (pd(:,:,time_stamp)) );
    assert( ~isreal(ph), 'complex-valued input required, ph' );

%     % Vector direction
%     [XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );
%     M = real( exp( 1i * angle(ph) ) ); 
%     N = imag( exp( 1i * angle(ph) ) );

    % 8 x 8 arrows
    [XX,YY] = meshgrid( 1:4:size(ph,2), 1:4:size(ph,1) );
    Y = conv2(ph,[1,1,1,1;1,1,1,1;1,1,1,1;1,1,1,1],'valid');
    Z = Y(1:4:end,1:4:end)/16
    M = real( exp( 1i * angle(Z) ) ); 
    N = imag( exp( 1i * angle(Z) ) );
     % quiver plot / Vector field
%     quiver( XX, YY, M, N, 0.7, 'k',  'linewidth', 1.5);

    % Streamlines Red
%   [XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );
    [XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );

    M = real( exp( 1i * angle(ph) ) ); 
    N = imag( exp( 1i * angle(ph) ) );
    [sx,sy] = meshgrid(4:16, 8:20); % Setup Steamline Origins in X: 1~32, Y: 1~32
%     [sx,sy] = meshgrid(1:32, 1:32); % Setup Steamline Origins in X: 1~32, Y: 1~32
	hlines_sens = streamline(stream2(XX,YY,M,N,sx,sy));
    set(hlines_sens,'LineWidth',1.5,'Color','r') % Red
    hold on;
    
    axis([-1 34 -1 34])

    
      drawnow 

%       % Capture each frame as an image 
%       frame = getframe(figure2); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       dtime = 1./frame_rate;
%       if i == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime); 
%       end 
        hold off;
    print('-dpng','-r300',txtTitle2)
    end
        close(figure2)