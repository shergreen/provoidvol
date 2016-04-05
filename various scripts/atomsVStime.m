%plots both import time and algorithm time as a function of Rprobe
%will look for a trendline

set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
figure('Color','white');

width = 8.37;     % Width in inches - adjust as necessary
height = 3.84;    % Height in inches - adjust as necessary

alw = 1.5;%0.75;    % AxesLineWidth 
fsz = 14;           % Fontsize 
lw = 1.5;           % LineWidth 
msz = 8;            % MarkerSize 


set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw); 

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','off'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
paperWidth = 0;
paperHeight = 1;
left = (defsize(1)- paperWidth)/2;
bottom = (defsize(2)- paperHeight)/2;
defsize = [left, bottom, 8, 10];
defsize = [0, 0, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

%the data is atoms, import time, algorithm time
hold on;
data = dlmread('FINALatomsVStime-Top500.txt');
%inertia = dlmread('inertiaFile-Top500.txt');
errorbar(data(:,1),data(:,3),data(:,4),'.')
%plot(inertia(:,1),inertia(:,2),'.')
%plot(inertia(:,1),inertia(:,3),'.')
%plot(inertia(:,1),inertia(:,4),'.')
xlabel('Number of atoms')
ylabel('Algorithm time (seconds)')
title('Linear algorithm time with protein atom count (0.5 Å lattice)')
hold off;

print('pub_quality_atoms_vs_time.png','-dpng','-r300');
