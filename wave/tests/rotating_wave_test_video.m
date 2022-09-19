% *WAVE*
%
% PHASE GRADIENT TEST SCRIPT     takes the phase gradient of a test
%                                   test signal (here, a 2d target wave)
%
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\test-data');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\analysis');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\plotting');

custFoldername = 'C:\Users\cybro\Dropbox\Matlab_codes\wave';

% clear all; clc %#ok<CLSCR>

% parameters
T = 1; %s
Fs = 1000; freq = 13.5; %Hz
image_size = 16; %px
pixel_spacing = 1; %a.u.
direction = +1; % +1/-1

% generate data
xf = zGrid;%generate_rotating_wave( image_size, 1/Fs, T, freq, direction );

% z-score data
xf = zscore_independent( xf );

% form analytic signal
xph = analytic_signal( xf );

% scale = 4;
% plot3((1:1001)-200, imag(reshape(xph(16,30,:), [1 1001])), real(reshape(xph(16,30,:), [1 1001])));
% hold on;
% plot3((1:1001)-200, scale*ones(1,1001), real(reshape(xph(16,30,:), [1 1001])));
% hold on;
% plot3((1:1001)-200, imag(reshape(xph(16,30,:), [1 1001])), -scale*ones(1,1001));
% axis([-200 800 -scale scale -scale scale ])
% xlabel('time (ms)');
% ylabel('Imaginary');
% zlabel('Real (z-score)');
% grid on;

% calculate instantaneous frequency 
[wt,signIF] = instantaneous_frequency( xph, Fs );

% calculate phase gradient
[pm,pd,dx,dy] = phase_gradient_complex_multiplication( xph, pixel_spacing, signIF );


% plot((1:1001)-200, reshape(pd(16,30,:), [1 1001]))
% hold on;
% plot((1:1001)-200, reshape(pd_unwrap(16,30,:), [1 1001]))
% 
% % plot(unwrap(reshape(pd(16,30,1:1001), [1 1001])))
% axis([-200 800 -3*pi 1.2*pi]);
% grid on;
% xlabel('time (ms)');
% ylabel('Degree (radian)');

   time_stamp = 220;
   dur = 0;
   time_ini = 205;
pd_unwrap = [];
for i=1:32
    for j=1:32
pd_unwrap(i,j,1:1001) = reshape(unwrap(reshape(pd(i,j,1:1001), [1 1001])), [1 1 1001]);
pd_unwrap(i,j,1:1001) = pd_unwrap(i,j,1:1001) - pd_unwrap(i,j,time_ini) + pd(i,j,time_ini);
    end
    i
end
% 
% for i = 1:32
% plot(reshape(pd_unwrap(i,3,:), [1 1001]));
% hold on;
% end
v_max = max(reshape(zGrid(:,:,1:1001), [32*32*1001 1]));
% v_max = max(reshape(pd_unwrap(:,:,1:1001), [32*32*1001 1]));
    v_min = -v_max;
    
    % Video output for 1s time segment from -0.2s to 0.8s
%Name of the video output file
v = VideoWriter(append(custFoldername, '\', txtTitle, '.avi'));
v.Quality = 95;
% Video frame rate & down sampling set up
v.FrameRate = 20; % How many frames per second for the video
open(v); 
size = 300;
ds = 1;
figure1=figure('Position', [150, 150, 150+size, 150+size]);

for  time_stamp = 150:2:551
  

   
%     figure1 = figure;
%         [C,h]= contourf(1:32,1:32,wt(:,:,time_stamp)*signIF,100);
% %       wt_HG = wt(:,:,time_stamp);
%         histogram(wt(:,:,time_stamp));
%         hold on;
%         histogram(wt_HG);
    [C,h]= contourf(1:32,1:32,zGrid(:,:,time_stamp),100);
%     [C,h]= contourf(1:32,1:32,pd_unwrap(:,:,time_stamp),100);
%     [C,h]= contourf(1:32,1:32,(pd_unwrap(:,:,time_stamp)-pd_unwrap(:,:,time_stamp-dur)),100);

    
    set(h,'LineColor','none')
    colormap('redblue');
    
      caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    
    txtTitle = append('Wave_Rat1_E4_whisker_beta', num2str(time_stamp-200), 'ms');
    title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([32 32 1]);
    set(gca,'YTick',0:4:31);
    set(gca,'YTickLabel',((0:4:31)-16)*0.15);
    set(gca,'XTick',0:4:31);
    set(gca,'XTickLabel',((0:4:31)-16)*0.15);
    set(gca,'FontSize',10);
    set(h,'FontSize',10);

    hold on;
% plot resulting vector field
if dur == 0
    plot_vector_field( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
else
    plot_vector_field( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
end

% Video mode
 pause(0.01);
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(v, frame);
%     % Save as image
%         h = gcf;
%         set(h,'PaperOrientation','landscape');
%         txt = append(custFoldername, '\', txtTitle, '.tif');
%         print('-dtiff','-r300',txt)
%         close(figure1)
        
end
        
hold off
close(v); % Saves the movie.
close(figure1);

%% VIDEO MODE
% % Video output for 1s time segment from -0.2s to 0.8s
% %Name of the video output file
% v = VideoWriter(append(custFoldername, '\', txtTitle, '.avi'));
% v.Quality = 95;
% Video frame rate & down sampling set up
% v.FrameRate = 20; % How many frames per second for the video
% open(v); 
% size = 300;
% ds = 1;
% figure1=figure('Position', [150, 150, 150+size, 150+size]);
% for i = 1:ds:fs+1
%     
%     pause(0.01);
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(v, frame);
%     end
%     
% hold off
% close(v); % Saves the movie.
% close(figure1);