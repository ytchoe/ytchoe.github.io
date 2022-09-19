% *WAVE*
%
% PHASE GRADIENT TEST SCRIPT     takes the phase gradient of a test
%                                   test signal (here, a 2d target wave)
%
txtName = 'Human_Feb28_2020_jump_correction_Pinky_vib'
custFoldername = 'C:\Users\cybro\Dropbox\Matlab_codes\Matlab_Feb28_OHSU\Plots_localization';
t_window = [100 200:5:300 400];

addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\test-data');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\analysis');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\wave\plotting');
addpath('C:\Users\cybro\Dropbox\Matlab_codes\export_fig-master');

% Color map
R = [zeros(1, 128) power(((0:127)/127),1)]
myColorMap = [zeros(256,1)+R', zeros(256,1), zeros(256,1)+flip(R)']

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

%    time_stamp = 220;
   dur = 0;
   time_ini = 205;
pd_unwrap = [];
for i=1:32
    for j=1:32
pd_unwrap(i,j,1:1000) = reshape(unwrap(reshape(pd(i,j,1:1000), [1 1000])), [1 1 1000]);
pd_unwrap(i,j,1:1000) = pd_unwrap(i,j,1:1000) - pd_unwrap(i,j,time_ini) + pd(i,j,time_ini);
    end
    i
end
% 
% for i = 1:32
% plot(reshape(pd_unwrap(i,3,:), [1 1001]));
% hold on;
% end

v_max = max(reshape(zGrid(:,:,201:300), [32*32*100 1]));
    v_min = -v_max;
% t_window = [50:50:200 210:10:300 350:50:500]
for  i = 1:length(t_window)
    time_stamp = t_window(i);
    figure1 = figure;
%         [C,h]= contourf(1:32,1:32,wt(:,:,time_stamp)*signIF,100);
% %       wt_HG = wt(:,:,time_stamp);
%         histogram(wt(:,:,time_stamp));
%         hold on;
%         histogram(wt_HG);
    colordef black;
    [C,h]= contourf(1:32,1:32,zGrid(:,:,time_stamp),100);

%     set(gca,'Color','k')

%     [C,h]= contourf(1:32,1:32,pd_unwrap(:,:,time_stamp),100); 
    set(h,'LineColor','none');
    colormap(myColorMap);
    axis([1 32 1 32]);
      caxis([v_min v_max]); 
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    
    txtTitle = append(txtName, '_', num2str(time_stamp-200), 'ms');
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    pbaspect([32 32 1]);
    set(gca,'YTick',0:4:31);
    set(gca,'YTickLabel',((0:4:31)-16)*0.15);
    set(gca,'XTick',0:4:31);
    set(gca,'XTickLabel',((0:4:31)-16)*0.15);
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
% if dur == 0
    plot_vector_field( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
%     hold on;
%     plot_vector_field2( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
%     hold on;
%     plot_vector_field4( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );

% else
%     plot_vector_field( exp( 1i .* (pd_unwrap(:,:,time_stamp)) ), 1 );
% end
    % Save as image
    set(gcf, 'Color', 'None');
        h = gcf;
        set(h,'PaperOrientation','landscape');
        txt = append(custFoldername, '\', txtTitle);
%          print('-dtiff','-r300',txt)
    export_fig(txt, '-png','-r300', '-transparent');
        
        close(figure1)
        
        end