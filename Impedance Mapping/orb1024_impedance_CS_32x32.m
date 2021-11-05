%%
%clear
%close all
% clc
% LGA1155_mapping_script

% File path for the CSV file
filepath = 'C:\JL_E3C\JL_E3B.csv';
% filepath = 'C:\16Feb2021_Awake\post-implant_pre-PhaseReversal_YT_B29S1L2_post-VPRO.csv';
% 
% imp_before = impedanceMag;
%imp_after = impedanceMag;

% Criteria for the good channels
ImpUpperBound = 500e3;
ImpLowerBound = 1e3;

[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');
impedanceMag = [];
impedancePhase = [];
%for i=1:1024
%impedanceMag(i) = amplifier_channels(i).electrode_impedance_magnitude;
%impedancePhase(i) = amplifier_channels(i).electrode_impedance_phase;
%end
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
impedancePhase = dlmread(filepath,',',[1,5,1024,5]);
%assuming that order of impedanceMag is A0-A63, B0-B63,...
%% load mapping files
load('impMagOrder');
load('eaglePadOrder');
load('padCoords');
load('elecCoords_CS_32x32');
load('ribbonCoords');
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
 elecX = elecCoords_CS_32x32(:,1);
 elecY = elecCoords_CS_32x32(:,2);
 ribX = ribbonCoords(:,1);
 ribY = ribbonCoords(:,2);
 
 counter = 0;
 lowImpedanceMag = [];
 lowImpedancePhase = [];
 lowPadX = [];
 lowPadY = [];
 lowElecX = [];
 lowElecY = [];
 k = [];
 for i = 1:length(impedanceMag)
     if (impedanceMag(i)> ImpUpperBound  || impedanceMag(i)< ImpLowerBound )
        %assume open or short
        colorVect(i,:) = [1 0 0];%red
     else
         k = [k i];
         counter = counter+1;
         lowImpedanceMag(counter) = impedanceMag(i);
         lowImpedancePhase(counter) = impedancePhase(i);
         lowPadX(counter) = padX(i);
         lowPadY(counter) = padY(i);
         lowElecX(counter) = elecX(i);
         lowElecY(counter) = elecY(i);
         colorVect(i,:) = [0 1 0];%green
     end
 end
disp(counter);
disp(mean(lowImpedanceMag));
sz = 10;

c = log10(impedanceMag);

avg = mean(lowImpedanceMag);
stdev = std(lowImpedanceMag);

p = lowImpedancePhase;
p_all = impedancePhase;
avg_p = mean(p);
stdev_p = std(p);
%%% Histogram of low impedances

f = figure('position', [100,100,1300,650]);

movegui('center');
subplot(2,3,2);

% c2 = [];
% m = 1;
% for i = k
% c2(m) = (imp_after(i)/imp_before(i));
% m = m + 1;
% end
% % edges2 = [0:0.01:2];
% histogram(c2);%, edges2);

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
h1 = histogram(p, edges);
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

%%% Mapping of all impedances on the Electrodes
subplot(2,3,4);
s = scatter(elecX(1:1024),elecY(1:1024),sz,log10(impedanceMag));
axis([-17 17 -17 17])
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

% Mapping of all impedances
%figure;
subplot(2,3,6);
s = scatter(padX(1:1024),padY(1:1024),sz,log10(impedanceMag));
axis([-21 21 -21 20])
box on
txt = ['Connector part, Count=' num2str(counter)];
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


% Mapping of ribbon cable metal leads impedances
%figure;
subplot(2,3,5);

X = log10(impedanceMag);
A = [ribX, X]
B = sortrows(A)
C = B(:,2);
C90 = rot90(C);
imagesc(C90)
%s = scatter(B(:,1),ribY,sz,C);
txt = ['Metal leads part, Count = ' num2str(counter)];
title(txt)
xlabel('x-position')
caxis(color_Range)
colorbar

h = gcf;
set(h,'PaperOrientation','landscape');
txt = [name '.pdf'];
saveas(gcf,txt)

%clearvars -except impedanceMag impedancePhase