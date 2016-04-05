%plots convergence of void/microvoid as a function of a_grid
d1aylH = csvread('1aylH-boundarywaterVSRprobe.txt');
d1cxyH = csvread('1cxyH-boundarywaterVSRprobe.txt');
d1mgtH = csvread('1mgtH-boundarywaterVSRprobe.txt');
d1moqH = csvread('1moqH-boundarywaterVSRprobe.txt');
d1notH = csvread('1notH-boundarywaterVSRprobe.txt');
d1qfmH = csvread('1qfmH-boundarywaterVSRprobe.txt');
d1qnfH = csvread('1qnfH-boundarywaterVSRprobe.txt');
d1qsaH = csvread('1qsaH-boundarywaterVSRprobe.txt');
d1ygeH = csvread('1ygeH-boundarywaterVSRprobe.txt');
d2driH = csvread('2driH-boundarywaterVSRprobe.txt');

data = zeros(16,3);

xlabel('R_{probe}')
ylabel('Total boundary solvent volume')
hold on;
for i=1:16
    data(i,1) = d1notH(i,1);
    data(i,2) = d1notH(i,2);
    data(i,3) = d1notH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1aylH(i,1);
    data(i,2) = d1aylH(i,2);
    data(i,3) = d1aylH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1cxyH(i,1);
    data(i,2) = d1cxyH(i,2);
    data(i,3) = d1cxyH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1mgtH(i,1);
    data(i,2) = d1mgtH(i,2);
    data(i,3) = d1mgtH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1moqH(i,1);
    data(i,2) = d1moqH(i,2);
    data(i,3) = d1moqH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1qfmH(i,1);
    data(i,2) = d1qfmH(i,2);
    data(i,3) = d1qfmH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1qnfH(i,1);
    data(i,2) = d1qnfH(i,2);
    data(i,3) = d1qnfH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1qsaH(i,1);
    data(i,2) = d1qsaH(i,2);
    data(i,3) = d1qsaH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d1ygeH(i,1);
    data(i,2) = d1ygeH(i,2);
    data(i,3) = d1ygeH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=1:16
    data(i,1) = d2driH(i,1);
    data(i,2) = d2driH(i,2);
    data(i,3) = d2driH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
%plotyy(d1notH(:,1),d1notH(:,2),d1notH(:,1),d1notH(:,4))
%errorbar(d1aylH(:,1),d1aylH(:,2),d1aylH(:,3))
%errorbar(d1cxyH(:,1),d1cxyH(:,2),d1cxyH(:,3))
%errorbar(d1mgtH(:,1),d1mgtH(:,2),d1mgtH(:,3))
%errorbar(d1moqH(:,1),d1moqH(:,2),d1moqH(:,3))
%errorbar(d1notH(:,1),d1notH(:,2),d1notH(:,3))
%errorbar(d1qfmH(:,1),d1qfmH(:,2),d1qfmH(:,3))
%errorbar(d1qnfH(:,1),d1qnfH(:,2),d1qnfH(:,3))
%errorbar(d1qsaH(:,1),d1qsaH(:,2),d1qsaH(:,3))
%errorbar(d1ygeH(:,1),d1ygeH(:,2),d1ygeH(:,3))
%errorbar(d2driH(:,1),d2driH(:,2),d2driH(:,3))


%plot(d1notH(:,1),d1notH(:,4))
hold off;