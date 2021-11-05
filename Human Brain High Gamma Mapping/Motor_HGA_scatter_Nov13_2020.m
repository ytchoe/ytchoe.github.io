    %% Plot HGA mapping
%     save('data_motor_HGA_Nov13_2020.mat', 'data_sigma');
    load('data_motor_HGA_Nov13_2020.mat');
    load('k_vector_CS_11132020.mat');
    load('elecCoords_CS_32x32.mat');
    elecX = elecCoords_CS_32x32(:,1);
    elecY = elecCoords_CS_32x32(:,2);

    f_size = 15;
	dtime = 1./20; 
    fs = 1e3; % Sampling rate
    c_max = 3; % Color range / Sigma high
    c_min = 1; % Color range / Sigma low


	filename = append('Motor_HGA_scatter_Nov13_2020.gif');;
    
    window_size_x = 800;
	window_size_y = window_size_x;
    pixel_size = 150;
    figure1=figure('Position', [0, 0, 0+window_size_x, 0+window_size_y]);
	set(gcf,'color','w');
    
    for i_gif = 850:5:(1600)
    
    scatter(elecX(k),elecY(k),pixel_size,data_sigma(k, i_gif),'filled', 'o'); 
    
    box on;
    ax=gca;
    ax.LineWidth=1;
    set(gca,'fontsize', f_size);
    set(gca,'fontweight', 'bold');
    caxis([c_min c_max]);

    h = colorbar;
	ylabel(h, 'HGA sigma (\sigma)');
    set(h,'FontSize',f_size)
    set(h,'FontSize',f_size)

% 	R = power((0:255)/255, 3);
%     myColorMap = [ones(256,1)-R', ones(256,1), ones(256,1)-R'];
%     colormap(gca, myColorMap);
    colormap('redblue');
    pbaspect([1 1 1]);
    axis([-17 17 -17 17]);
    xlabel('mm');
    ylabel('mm');

    txt = [num2str(i_gif-1000, '%d ms')];
    title(txt);        

     % Capture the plot as an image 
      frame = getframe(figure1); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256);
      
    % Write to the GIF File 
      if i_gif == 850 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',dtime); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',dtime); 
      end 
        hold off;
    
    end
    close(figure1);