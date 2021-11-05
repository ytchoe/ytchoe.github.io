%% Impedance mapping function
% Last modified by Youngbin Tchoe, 8/23/2021
function k = orb1024_impedance_ratHD_yb(filepath, ImpLowerBound, ImpUpperBound)

[flocation,name,ext] = fileparts(filepath)
fid = fopen(filepath, 'rt');

% Reading impedance magnitude and phase @ 1 kHz
impedanceMag = [];
impedancePhase = [];
impedanceMag = dlmread(filepath,',',[1,4,1024,4]);
impedancePhase = dlmread(filepath,',',[1,5,1024,5]);

%% load mapping files
load('impMagOrder');
load('eaglePadOrder');
load('padCoords');
load('elecCoords_ratHD');
load('ribbonCoords');
%% Load coordinate information
 colorVect = zeros(1027,3);
 padX = newCoords(:,1);
 padY = newCoords(:,2);
 elecX = elecCoords_ratHD(:,1);
 elecY = elecCoords_ratHD(:,2);
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

edges = [2:0.05:7];
h1 = histogram(c, edges);
h1.FaceColor = [1,0,0];
xticks([2 3 4 5 6 7])
xticklabels({'100','1K','10K','100K', '1M', '10M'})
txt = ['Impedance Magnitude Histogram'];
title(txt);
xlabel('Impedance Magnitude @ 1 kHz')
ylabel('Counts')

%%% Phase: Histogram of low impedances
subplot(2,3,3);
edges = [-180:3:180];
h1 = histogram(p, edges);
h1.FaceColor = [0,1,0];

title('Phase Histogram')
xticks([-180 -135 -90 -45 0 45 90 135 180])
xtickangle(45)
xlabel('Phase @ 1 kHz')
ylabel('Counts')


%%% Mapping of all impedances on the Electrodes
subplot(2,3,4);
s = scatter(elecX(1:1024),elecY(1:1024),sz,log10(impedanceMag));
axis([-2.7 2.7 -2.7 2.7])
box on
txt = ['Electrodes part, Count =' num2str(counter)];
title(txt, 'Interpreter', 'none');
xlabel('x-position')
ylabel('y-position')
s.Marker = 'o';
color_Range = [3.5 6.5];
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

s.LineWidth = 0.75;
caxis(color_Range)
colorbar


h3 = subplot(2,3,1);
txt = [name char(10) char(10) 'Working channel (<' num2str(int32(ImpUpperBound/1000)) 'kOhm & >' num2str(int32(ImpLowerBound/1000)) 'kOhm @ 1kHz)' char(10) ' = ' num2str(counter) char(10) 'Impedance average = ' num2str(avg, 2) ' Ohm' char(10) 'Standard deviation = ' num2str(stdev, 2) ' Ohm'];
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

%% 'k' saves the channels in reliable impedance magnitude range and excludes
k = [];
j = 1;
 for i = 1:1024
     
     if(impedanceMag(i)<ImpUpperBound && impedanceMag(i)>ImpLowerBound) %&& (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1))&& (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1))&& (~(abs(elecX(i)+1.425) < 0.1 && abs(elecY(i)+2.325) < 0.1)) && (~(abs(elecX(i)-0.075) < 0.1 && abs(elecY(i)+1.725) < 0.1)) && (~(abs(elecX(i)+0.675) < 0.1 && abs(elecY(i)+2.025) < 0.1)) && (~((abs(elecX(i)-0.15) < 0.15) && (elecY(i)>2.3))) && (~(abs(elecX(i)-0.525) < 0.1 && abs(elecY(i)-2.325) < 0.1)) && (~(abs(elecX(i)-2.325) < 0.1 && abs(elecY(i)-2.175) < 0.1)) && (~(abs(elecX(i)-2.175) < 0.1 && abs(elecY(i)-2.175) < 0.1))
         k(1,j) = i;
         j = j+1;
     end
 end
 length(k)
 
 return