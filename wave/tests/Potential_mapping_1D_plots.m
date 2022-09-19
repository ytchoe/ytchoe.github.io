load('C:\Plots_Rat2_RMS\TaCMR_Electrical_stim__hindlimb_2p55mA_200317_160038.mat');
load('elecCoords_ratHD');
elecX = elecCoords_ratHD(:,1);
elecY = elecCoords_ratHD(:,2);
load('k_vector_rat2.mat');
% load('High_Gamma_low_order.mat');

custFoldername = 'C:\Users\cybro\Dropbox\Matlab_codes';
% generate data
data = ampDataTaCMR;
fs = 20000;

time_stamp = 220;

% v_max = max(reshape(data(k,4011:5010), [length(k)*1000 1]))
% v_max = max(max(data(k,4401:4441)))
v_max = 50;
v_min = -v_max;
% Color map
R = [zeros(1, 128) (0:127)/127];
myColorMap = [R', zeros(256,1), flip(R')];

for  time_stamp = 8000:20:8000
    figure1 = figure;
%         [C,h]= contourf(1:32,1:32,wt(:,:,time_stamp)*signIF,100);
% %       wt_HG = wt(:,:,time_stamp);
%         histogram(wt(:,:,time_stamp));
%         hold on;
%         histogram(wt_HG);
    colordef black
    scatter(elecX(k),elecY(k),120,data(k,time_stamp),'filled', 's');
    box on;
    colormap(myColorMap);
    set(gca,'Color','w')
    caxis([v_min v_max]); 
    axis([-2.5 2.5 -2.5 2.5]);
    pbaspect([32 32 1]);
    h = colorbar;
    xlabel('mm');
    ylabel('mm');
    txtTitle = append('Rat2_Hindlimb_Electric_Stim_2.55mA_', num2str((time_stamp-4000)/fs*1000, '%4.0f'), 'ms');
    ht = title(txtTitle,'interpreter', 'none');
    ylabel(h, 'Voltage (\muV)');
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    set(gca,'color','none')
    set(gca,'linewidth',1)
    set(h,'linewidth',1)
    set(h,'FontSize',10);
    set(h,'FontWeight','bold');
    set(ht,'FontWeight','bold');


        % Save as image
        set(gcf, 'Color', 'None')
        h = gcf;
        set(h,'PaperOrientation','landscape');
        txt = append(custFoldername, '\', txtTitle);
%       print('-dtiff','-r300',txt)
        export_fig(txt, '-png','-r300', '-transparent');
        close(figure1)
        
end

%% Plot along y = Const.
    vMargin_up = 50;
    vMargin_down = 100;
    vOffset = 20;
    name = 'Rat2_Electric_Stim_2.55mA';
    for scan = 1:32
    n_y = [];
    j = 1;
    y_pos = 17*0.15-0.075 -scan*0.15
    for i = 1:length(k)
        if (abs(elecY(k(i)) - y_pos)) < 0.1
            n_y(1, j) = elecX(k(i));
            n_y(2, j) = k(i);
            j = j + 1;
        end
    end
    j
%     Sorting n_y
    n_y = transpose(sortrows(transpose(n_y), 'descend')); %'descend'
    
%   Plot along the line in n_y datas
    time1 = -0.1
    time2 = 0.2
    range = (time1+0.2)*fs:(time2+0.2)*fs;
    figure2 = figure('Position', [50, 50, 150+500, 150+800]);
    for i = 1:length(n_y(2,:))
    plot((range)/fs-0.2, vOffset*i+data(n_y(2,i),range),'LineWidth',2);
    line([0 0], ylim, 'color', 'white')
    xlabel('Time (s)');
    txt2 = append(name, '_y_', num2str(elecY(n_y(2, 1))));
    h = title(txt2,'interpreter', 'none');
    ylabel('Voltage (uV)');
    axis([time1 time2 -vMargin_down vMargin_up+vOffset*32]);%length(n_y(2,:))]);
    set(gca,'FontSize',15)
    set(gca,'FontWeight','bold');
    set(h,'FontWeight','bold');
    set(gca,'color','none')
    set(gca,'linewidth',1)
    hold on;
    end
%     Save as image
    h = gcf;
    set(h,'PaperOrientation','landscape');
    txt = append('Waveform_', name, '_y_', num2str(elecY(n_y(2, 1))));
%     saveas(gcf,txt)
    export_fig(txt, '-png','-r300', '-transparent');
    close(figure2);
    end

%% Plot along x = Const.
    vMargin_up = 50;
    vMargin_down = 100;
    vOffset = 20;
    name = 'Rat2_Electric_Stim_2.55mA';
    for scan = 1:32
    n_x = [];
    j = 1;
    x_pos = 17*0.15-0.075 -scan*0.15
    for i = 1:length(k)
        if (abs(elecX(k(i)) - x_pos)) < 0.1
            n_x(1, j) = elecY(k(i));
            n_x(2, j) = k(i);
            j = j + 1;
        end
    end
    j
%     Sorting n_x
    n_x = transpose(sortrows(transpose(n_x), 'descend')); %'descend'
    
%   Plot along the line in n_y datas
    time1 = -0.1
    time2 = 0.2
    range = (time1+0.2)*fs:(time2+0.2)*fs;
    figure2 = figure('Position', [50, 50, 150+500, 150+800]);
    for i = 1:length(n_x(2,:))
    plot(((range)/fs-0.2)*1000, vOffset*(33-i)+data(n_x(2,i),range),'LineWidth',2);
    line([0 0], ylim, 'color', 'white')
    xlabel('Time (ms)');
    txt2 = append(name, '_x_', num2str(elecX(n_x(2, 1))));
    h = title(txt2,'interpreter', 'none');
    ylabel('Voltage (uV)');
    axis([-50 100 -vMargin_down vMargin_up+vOffset*32]);%length(n_y(2,:))]);
    set(gca,'FontSize',15)
    set(gca,'FontWeight','bold');
    set(h,'FontWeight','bold');
    set(gca,'color','none')
    set(gca,'linewidth',1)
    hold on;
    end
%     Save as image
    h = gcf;
    set(h,'PaperOrientation','landscape');
    txt = append('Waveform_', name, '_x_', num2str(elecX(n_x(2, 1))));
%     saveas(gcf,txt)
    export_fig(txt, '-png','-r300', '-transparent');
    close(figure2);
end