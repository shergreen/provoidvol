%plots convergence of void/microvoid as a function of a_grid
proteins = dir('*.txt');
xlabel('R_{probe}')
ylabel('Total void volume (angstroms cubed)')
hold on;
for i = 1:length(proteins)
    data = csvread(proteins(i).name);
    errorbar(data(:,1),data(:,2),data(:,3))
end
hold off;