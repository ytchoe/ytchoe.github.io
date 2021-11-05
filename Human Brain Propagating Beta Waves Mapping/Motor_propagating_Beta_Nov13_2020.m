    frame_rate = 10; % Frame rate of GIF
    dtime = 1/frame_rate; % (s) Delay time between GIF frames. 1/dtime: 

    load('prop_beta_wave_13Nov2020.mat'); % data_beta: Beta wave potential / data_pd: Phase gradient direction
    load('data_motor_HGA_Nov13_2020.mat'); % HGA sigma
    load('k_vector_CS_11132020.mat');
    load('elecCoords_CS_32x32.mat');
    elecX = elecCoords_CS_32x32(:,1);
    elecY = elecCoords_CS_32x32(:,2);

    f_size = 15; % Font size
    fs = 1e3; % Sampling rate
    c_max = 100; % Color range / Sigma high
    c_min = -100; % Color range / Sigma low

	filename = append('Motor_beta_after_Nov13_2020.gif');;
    
    % Time window
	time_zero = 1000; % t = 0, motion starts
    time_step = 1; % ms
    time_first = 350; % ms, prior to motion
    time_last = 450; % ms, after the motion
    t_window = [time_zero+time_first:time_step:time_zero+time_last];
    
    pixel_size = 150;
    x_window_size = 800; % Window size
    y_window_size = 800;
    figure1 = figure('Position', [0, 0, 0+x_window_size, 0+y_window_size]);
    set(gcf,'color','w');
    
for  i = 1:length(t_window)
    time_stamp = t_window(i);
    
    % Background beta map
	[C,h]= contourf(1:32,1:32,data_beta(:,:,time_stamp),100);
    set(h,'LineColor','none');
    
    pbaspect([32 32 1]);
    set(gca,'YTick',0:4:31);
    set(gca,'YTickLabel',((0:4:31)-16));
    set(gca,'XTick',0:4:31);
    set(gca,'XTickLabel',((0:4:31)-16));
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    set(gca,'color','none');
    set(gca,'linewidth',1);

    
    box on;
    ax=gca;
    ax.LineWidth=1;
    set(gca,'fontsize', f_size);
    set(gca,'fontweight', 'bold');
    caxis([c_min c_max]);

    hc = colorbar;
    ylabel(hc, 'Voltage (\muV)');
    set(hc,'FontSize',f_size)
	set(hc,'linewidth',1);
    set(hc,'FontWeight','bold');
%     set(ht,'FontWeight','bold');
%     set(h,'FontSize',f_size)

    % Color map
% 	R = power((0:255)/255, 3);
% 	myColorMap = [ones(256,1), ones(256,1)-R', ones(256,1)-R']; % Red colormap
%     myColorMap = [ones(256,1)-R', ones(256,1), ones(256,1)-R']; % Green colormap
%     colormap(gca, myColorMap);
    colormap('redblue');
    pbaspect([1 1 1]);
    xlabel('mm');
    ylabel('mm');

    % Plot Title
    t = time_stamp-time_zero;
    if ((t > -500) && (t < -80)) txt_stage = 'Preparing for movement';
    elseif (t > -80) && (t < 190) txt_stage = 'Hand grabbing motion';
    elseif (t > 190) && (t < 290) txt_stage = 'Completing the movement';
	elseif (t > 290) && (t < 500) txt_stage = 'Completed the movement';
    else txt_stage = ' ';
    end
    
    txtTitle = append( {txt_stage, append(num2str(t, '%d'), 'ms')} );    
    ht = title(txtTitle,'interpreter', 'none');
	set(ht,'FontWeight','bold');

    hold on;
    
    % Plot Vector Field
    ph = exp( 1i .* (data_pd(:,:,time_stamp)) );
    assert( ~isreal(ph), 'complex-valued input required, ph' );

    % 
    [XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );
    M = real( exp( 1i * angle(ph) ) ); 
    N = imag( exp( 1i * angle(ph) ) );
    
%     % Streamlines Red
%     [sx,sy] = meshgrid(6:24, 3:10); % Setup Steamline Origins in X: 1~32, Y: 1~32
% 	hlines_sens = streamline(stream2(XX,YY,M,N,sx, sy));
%     set(hlines_sens,'LineWidth',2,'Color','r') % Red
%     hold on;
%     
%     % Streamlines Blue
%     [sx,sy] = meshgrid(16:32, 16:24); % Setup Steamline Origins in X: 1~32, Y: 1~32
% 	hlines_motor = streamline(stream2(XX,YY,M,N,sx,sy));
%     set(hlines_motor,'LineWidth',2,'Color','b') % Blue
%     hold on;
%     
    % quiver plot
    quiver( XX, YY, M, N, 0.75, 'k',  'linewidth', 1);
%     axis equal
    axis([1 32 1 32]);

    drawnow 
    
      % GIF part
      % Capture the plot as an image 
      frame = getframe(figure1); 
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
        close(figure1)