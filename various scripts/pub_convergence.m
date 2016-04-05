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

%plots convergence of void/microvoid as a function of a_grid
hold on;
proteins = dir('*-gridCheck.txt');
for i=1:length(proteins)
raw_data = dlmread(proteins(i).name);
baseline = raw_data(1,2); %this is the total cavity for a metric
errorbar(raw_data(:,1),raw_data(:,2)/baseline,raw_data(:,3)/baseline);
end
xlabel('lattice constant (‎Å)')
ylabel('normalized total cavity volume')
title('Convergence with decreasing lattice constant')

hold off;

print('pub_quality_convergence.png','-dpng','-r300');
