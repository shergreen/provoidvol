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
proteins = dir('*-gridCheck.txt');
for i=1:length(proteins)
raw_data = dlmread(proteins(i).name);
baseline = raw_data(1,4); %this is the maximum length of time for a metric
inverse_cubed = (1 ./raw_data(:,1)).^3;
plot(inverse_cubed,raw_data(:,4),'-x');
end
xlabel('1/lattice constant^3 (‎Å^{-3})')
ylabel('Normalized total algorithm time')
title('Algorithm time with (1/lattice)^3')

hold off;

print('pub_quality_lattice_time.png','-dpng','-r300');
