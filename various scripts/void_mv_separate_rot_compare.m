%this plots the error over orientations separately for microvoid and cavity
%for the top 10

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

%protein size(x) avg void (std dev) avg mv (std dev)

data = dlmread('data_file.txt');
hold on;
subplot(2,1,1);
title('Normalized cavity and microvoid error over 100 orientations');
errorbar(data(:,1),data(:,3)./data(:,3),data(:,4)./data(:,3),'b+');
ylabel('cavity');
xlabel('Protein  volume (Å^3)');
subplot(2,1,2);
errorbar(data(:,1),data(:,5)./data(:,5),data(:,6)./data(:,5),'r+');
ylabel('microvoid');
hold off;



print('rotation_error_mv_and_void.png','-dpng','-r300');