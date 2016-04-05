%compares the vdw distribution to random rotation stats
vdw_data = dlmread('Bondi_vdw_error_nonames.txt');
rot_data = dlmread('RotationIndependence.txt');
volumes = dlmread('top10_vols.txt');

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


hold on;
errorbar(volumes,vdw_data(:,1)./vdw_data(:,1),vdw_data(:,2)./vdw_data(:,1),'b+');
errorbar(volumes,rot_data(:,1)./rot_data(:,1),rot_data(:,2)./rot_data(:,1),'r+');
title('Normalized error from vdW perturbations vs. random rotations');
xlabel('Protein volume (Å^3)');
ylabel('Normalized total cavity + microvoid');
legend('vdW perturbations','Random rotations');
hold off;

print('vdwVSrotations.png','-dpng','-r300');
