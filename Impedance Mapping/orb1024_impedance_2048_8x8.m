%% LGA1155_mapping_script
% Open figure
f = figure('position', [100,100,1300,650]);

%% load mapping files
load('impMagOrder');
load('eaglePadOrder');
load('padCoords');
load('elecCoords_8x8_L');
load('elecCoords_8x8_R');
load('ribbonCoords');

% Criteria for the good channels
ImpUpperBound = 200e3;
ImpLowerBound = 1e3;

% File path for the CSV file
% filepath_L = 'H:\Human\OHSU_awake_motor_Sep_08_2020\08Sept2020_Awake\Implantation1_L.csv';
% filepath_R = 'H:\Human\OHSU_awake_motor_Sep_08_2020\08Sept2020_Awake\Implantation1_R.csv';
filepath_L = 'H:\Human\OHSU_awake_motor_Sep_08_2020\08Sept2020_Awake\Post_implant2_L.csv';
filepath_R = 'H:\Human\OHSU_awake_motor_Sep_08_2020\08Sept2020_Awake\Post_implnat2_R.csv';


[flocation,name,ext] = fileparts(filepath_L)
fid = fopen(filepath_L, 'rt');
impedanceMag_L = [];
impedancePhase_L = [];

[flocation,name,ext] = fileparts(filepath_R)
fid = fopen(filepath_R, 'rt');
impedanceMag_R = [];
impedancePhase_R = [];
%for i=1:1024
%impedanceMag(i) = amplifier_channels(i).electrode_impedance_magnitude;
%impedancePhase(i) = amplifier_channels(i).electrode_impedance_phase;
%end
impedanceMag_L = dlmread(filepath_L,',',[1,4,1024,4]);
impedancePhase_L = dlmread(filepath_L,',',[1,5,1024,5]);
impedanceMag_R = dlmread(filepath_R,',',[1,4,1024,4]);
impedancePhase_R = dlmread(filepath_R,',',[1,5,1024,5]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...

%%
%padName matches coords, and this is what we are plotting. So we need to
%get impedanceMag to match padName's order
%Basically, signalName and padMatch indicies match up, so we need to use
%this to change impedanceMag's order to match padName

%first, we change impedanceMag's order to match padMatch (or signalName).
%Then we change the order to match padName
% tmpImpMag = zeros(1,length(signalName)-3);
% 
% for i = 1:(length(signalName)-3)
%     curName = impMagOrder{i};
%     curIndex = find(endsWith(signalName(:),curName));
%     tmpImpMag(i) = impedanceMag(curIndex);
% end
% 
% %now we need to map to the padName order from padMatch order (which is what
% %we're currently in)
% for i = 1:(length(padName)-3)
%     curName = padName{i};
%     curIndex = find(contains(padMatch{1},curName));
%     impedanceMag(i) = tmpImpMag(i);
% end

 colorVect = zeros(1027,3);
 padX = newCoords(:,1);
 padY = newCoords(:,2);
 elecX_L = elecCoords_8x8_L(:,1);
 elecY_L = elecCoords_8x8_L(:,2);
 elecX_R = elecCoords_8x8_R(:,1);
 elecY_R = elecCoords_8x8_R(:,2);
 ribX = ribbonCoords(:,1);
 ribY = ribbonCoords(:,2);
 
 counter_L = 0;
 lowImpedanceMag_L = [];
 lowImpedancePhase_L = [];
%  lowPadX = [];
%  lowPadY = [];
%  lowElecX = [];
%  lowElecY = [];
 for i = 1:length(impedanceMag_L)
     if (impedanceMag_L(i)> ImpUpperBound  || impedanceMag_L(i)< ImpLowerBound )
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter_L = counter_L+1;
         lowImpedanceMag_L(counter_L) = impedanceMag_L(i);
         lowImpedancePhase_L(counter_L) = impedancePhase_L(i);
%          lowPadX(counter_L) = padX(i);
%          lowPadY(counter_L) = padY(i);
%          lowElecX(counter_L) = elecX_L(i);
%          lowElecY(counter_L) = elecY_L(i);
         colorVect(i,:) = [0 1 0];%green
     end
 end
 
  counter_R = 0;
 lowImpedanceMag_R = [];
 lowImpedancePhase_R = [];
%  lowPadX = [];
%  lowPadY = [];
%  lowElecX = [];
%  lowElecY = [];
 for i = 1:length(impedanceMag_R)
     if (impedanceMag_R(i)> ImpUpperBound  || impedanceMag_R(i)< ImpLowerBound )
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         counter_R = counter_R+1;
         lowImpedanceMag_R(counter_R) = impedanceMag_R(i);
         lowImpedancePhase_R(counter_R) = impedancePhase_R(i);
%          lowPadX(counter_L) = padX(i);
%          lowPadY(counter_L) = padY(i);
%          lowElecX(counter_L) = elecX_L(i);
%          lowElecY(counter_L) = elecY_L(i);
         colorVect(i,:) = [0 1 0];%green
     end
 end
disp(counter_R);
disp(mean(lowImpedanceMag_R));
sz = 10;

avg = mean([lowImpedanceMag_L lowImpedanceMag_R]);
stdev = std([lowImpedanceMag_L lowImpedanceMag_R]);

counter = counter_L + counter_R;

p = lowImpedancePhase_L;
p_all = [impedancePhase_L impedancePhase_R];
avg_p = mean(p);
stdev_p = std(p);
%% Histogram of low impedances
c = log10([impedanceMag_L impedanceMag_R]);


% movegui('center');
subplot(2,3,2);

edges = [2:0.05:7];
h1 = histogram(c, edges);
h1.FaceColor = [1,0,0];
xticks([2 3 4 5 6 7])
xticklabels({'100','1K','10K','100K', '1M', '10M'})
txt = ['Impedance Magnitude Histogram'];
title(txt);
xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')

%dim = [0.13 0.8 0.3 0.1];
%annotation('textbox',dim,'String',txt,'FitBoxToText','on');

%%% Phase: Histogram of low impedances
%figure;
subplot(2,3,3);
edges = [-180:3:180];
h1 = histogram(p_all, edges);
h1.FaceColor = [0,1,0];
%hold on
%h2 = histogram(p_all, edges);
%h2.FaceColor = [0,0,1];
title('Phase Histogram')
xticks([-180 -135 -90 -45 0 45 90 135 180])
xtickangle(45)
xlabel('Phase @ 1 kHz')
ylabel('Counts')
%%% Histogram of all impedances
%figure;
%edges = [0:0.25:8];
%histogram(log10(impedanceMag), edges);

%% Mapping of all impedances on the Electrodes
subplot(2,3,4);
s = scatter(elecX_L(1:1024),elecY_L(1:1024),sz,log10(impedanceMag_L));
hold on;
s = scatter(elecX_R(1:1024),elecY_R(1:1024),sz,log10(impedanceMag_R));
axis([-42 42 -42 42])
box on
txt = ['Electrodes part, Count =' num2str(counter)];
title(txt, 'Interpreter', 'none');
xlabel('x-position')
ylabel('y-position')
s.Marker = 'o';
color_Range = [3.5 6.5];
%s.MarkerEdgeColor = 'flat';
%s.MarkerFaceColor = 'flat';
s.LineWidth = 0.75;
caxis(color_Range)
colorbar

%% Mapping of impedance on the connector region_L
%figure;
subplot(2,3,5);
s = scatter(padX(1:1024),padY(1:1024),sz,log10(impedanceMag_L));
axis([-21 21 -21 20])
box on
txt = ['Left Connector, Count=' num2str(counter_L)];
title(txt)
xlabel('x-position')
ylabel('y-position')
s.Marker = 'o';
%s.MarkerEdgeColor = 'flat';
%s.MarkerFaceColor = 'flat';
s.LineWidth = 0.75;
caxis(color_Range)
colorbar
%colormap([0 1 0;0 1 0;0 1 0;0 1 0;1 0 0;1 0 0;1 0 0]);

% scatter3(padX(1:1024),padY(1:1024),impedanceMag)
% set(gca,'zscale','log')

%% Mapping of impedance on the connector region_L
%figure;
subplot(2,3,6);
s = scatter(padX(1:1024),padY(1:1024),sz,log10(impedanceMag_R));
axis([-21 21 -21 20])
box on
txt = ['Right Connector, Count=' num2str(counter_R)];
title(txt)
xlabel('x-position')
ylabel('y-position')
s.Marker = 'o';
%s.MarkerEdgeColor = 'flat';
%s.MarkerFaceColor = 'flat';
s.LineWidth = 0.75;
caxis(color_Range)
colorbar
%colormap([0 1 0;0 1 0;0 1 0;0 1 0;1 0 0;1 0 0;1 0 0]);

% scatter3(padX(1:1024),padY(1:1024),impedanceMag)
% set(gca,'zscale','log')

%%% blank
%figure;
h3 = subplot(2,3,1);
txt = [name char(10) char(10) 'Working channel (<' num2str(int32(ImpUpperBound/1000)) 'kOhm & >' num2str(int32(ImpLowerBound/1000)) 'kOhm @ 1kHz)' char(10) ' = ' num2str(counter) char(10) 'Impedance average = ' num2str(avg, 2) ' Ohm' char(10) 'Standard deviation = ' num2str(stdev, 2) ' Ohm'];
%title(txt, 'Interpreter', 'none');
t = uicontrol(f,'Style','text','String','Select a data set.');
t.String = txt;
t.Position = [120 400 350 200];
t.FontSize = 12;
t.BackgroundColor = 'white'

set(h3,'Visible','off');


% % Mapping of ribbon cable metal leads impedances
% %figure;
% subplot(2,3,5);
% 
% X = log10(impedanceMag_L);
% A = [ribX, X]
% B = sortrows(A)
% C = B(:,2);
% C90 = rot90(C);
% imagesc(C90)
% %s = scatter(B(:,1),ribY,sz,C);
% txt = ['Metal leads part, Count = ' num2str(counter)];
% title(txt)
% xlabel('x-position')
% caxis(color_Range)
% colorbar
% 
h = gcf;
set(h,'PaperOrientation','landscape');
txt = [name '.pdf'];
saveas(gcf,txt)
% 


%% 'k' saves the channels in reliable impedance magnitude range and excludes
% the data with too high high gamma values at a selected time segment (which looked like a bad
% % channel)
 % Left side

k_L = [];
j = 1;
 for i = 1:1024
     if(impedanceMag_L(i)<ImpUpperBound && impedanceMag_L(i)>ImpLowerBound)
         k_L(1,j) = i;
         j = j+1;
     end
 end
 % Right side
k_R = [];
j = 1;
 for i = 1:1024
     if(impedanceMag_R(i)<ImpUpperBound && impedanceMag_R(i)>ImpLowerBound)
         k_R(1,j) = i;
         j = j+1;
     end
 end
%  
%  save('k_L_vector_Sep08_2020_1st.mat', 'k_L')
%  save('k_R_vector_Sep08_2020_1st.mat', 'k_R')
%  
%  clearvars -except impedanceMag_L impedanceMag_R k_L k_R
%  