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


%print('filename','-dpng','-r300');

a_grid = 0.5;
proteins = dir('*-microvoidHist.txt');
names = dir('*-totalVolume.txt');

% hold on;
% for i = 1:length(proteins)
%     data = csvread(proteins(i).name);
%     totalVolume = csvread(names(i).name);
%     data(:,1) = data(:,1) .* (0.5*0.5*0.5) ./ (totalVolume.^(1/3));
%     plot(log10(data(:,1)),log10(data(:,2)),'.')
% end
% hold off;

%%next we want to get the total volumes so we can scale and see if they
%%converge to one

%%now we want to average all the cluster sizes with one frequency

%we also want to generate a line that fits the entire data set
totalData = [];

hold on;
for i = 1:length(proteins)
    data = csvread(proteins(i).name);
    
    averagedData = zeros(size(data,1),size(data,2));
    for j = 1:length(data)
        clusterAmount = data(j,2);
        count = 0;
        sum = 0;
        if any(abs(averagedData(:,2)-clusterAmount)<1e-10) %if we have already dealt with this cluster quantity
            %don't do anything, we've done it already
        elseif (clusterAmount < 10) %if the cluster amount is less than our threshold
            %we dont do anything, we dont want this data
        else
            %take average over all with this quantity
            for k=1:length(data)
                if (abs(data(k,2)-clusterAmount)<1e-10) %if the cluster quantity is the same as one in question
                    sum = sum + data(k,1);
                    count = count + 1;
                end
            end
            averagedData(j,1) = sum / count; %average cluster size for this quantity
            averagedData(j,2) = clusterAmount;
        end
    end
    totalVolume = csvread(names(i).name);
    data(:,1) = data(:,1) .* (a_grid*a_grid*a_grid) ./ (totalVolume.^(1/3));
    averagedData(:,1) = averagedData(:,1) ./ (totalVolume.^(1/3));
    %totalData = vertcat(totalData,averagedData);
    %plot(log10(data(:,1)),log10(data(:,2)),'.')
    plot(log10(averagedData(:,1)),log10(averagedData(:,2))-2,'.')
end
%totalData(all(totalData==0,2),:)=[];
%totalData(:,2) = totalData(:,2) * a_grid;
%plot(log10(totalData(:,1)),log10(totalData(:,2))-2,'.')
hold off;