% 	load('data_stim_evoked_4th_epilepsy_Sep25_2020.mat');
% 	load('data2_stim_evoked_4th_epilepsy_Sep25_2020.mat');
% 	save('data2_stim_evoked_4th_epilepsy_Sep25_2020.mat', 'data_scatter_amp', 'data_amp');
    load('data_spontaneous_epilepsy_Sep25_2020.mat');
    
    load('k_vector_Sep25_2020.mat');
    load('elecCoords_CS_32x32.mat');
    elecX = elecCoords_CS_32x32(:,1);
    elecY = elecCoords_CS_32x32(:,2);
%         save('data_spontaneous_epilepsy_Sep25_2020.mat', 'data', 'data_amp', 'data_scatter_amp', 'data_pd');

    frame_rate = 24; % Frame rate of GIF
    dtime = 1/frame_rate; % (s) Delay time between GIF frames. 1/dtime: 
    fName = 'Amplitude2_scatter_25Sep2020_Spontaneous_v2'; % GIF filename
    filename = append(fName, '.gif');
    v_max = 100; % Potential range max in uV
    v_min = 0; % Potential range min in uV
    pixel_size = 150; % Sactter plot dot size
    f_size = 15; % Font size
    
    % Window size
    x_window_size = 800; 
    y_window_size = 800;
    figure1 = figure('Position', [0, 0, 0+x_window_size, 0+y_window_size]);
    set(gcf,'color','w');
    
    % Time window
	time_zero = 0; % t = 0, motion starts
    time_step = 10; % ms
    time_first = 1-time_zero; % ms
    time_last = 4001-time_zero; % ms
    t_window = [time_zero+time_first:time_step:time_zero+time_last];
%     
%     data_scatter_amp = data_amp;
    plot(data_scatter_amp(1,:))
for  i = 1:length(t_window)
    time_stamp = t_window(i);
    
    % Scatter plot
	scatter(elecX(k),elecY(k), pixel_size, data_scatter_amp(k, time_stamp),'filled', 'o'); 
    colormap(flipud(hot));
%     colormap(jet);
	lim = 17;
    axis([-lim lim -lim lim]);
    caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    
    txtTitle = append( num2str(time_stamp-time_zero), ' ms' );
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