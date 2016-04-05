%plots convergence of void/microvoid as a function of a_grid
d1aylH = csvread('1aylH-Rprobe1.0-gridCheck.txt');
d1cxyH = csvread('1cxyH-Rprobe1.0-gridCheck.txt');
d1mgtH = csvread('1mgtH-Rprobe1.0-gridCheck.txt');
d1moqH = csvread('1moqH-Rprobe1.0-gridCheck.txt');
d1notH = csvread('1notH-Rprobe1.0-gridCheck.txt');
d1qfmH = csvread('1qfmH-Rprobe1.0-gridCheck.txt');
d1qnfH = csvread('1qnfH-Rprobe1.0-gridCheck.txt');
d1qsaH = csvread('1qsaH-Rprobe1.0-gridCheck.txt');
d1ygeH = csvread('1ygeH-Rprobe1.0-gridCheck.txt');
d2driH = csvread('2driH-Rprobe1.0-gridCheck.txt');

data = zeros(18,3);
k = 100;


xlabel('a_{grid}')
ylabel('Change in void between a_{grid} sizes')
hold on;
for i=18:-1:1
    data(i,1) = d1notH(i+1,1);
    data(i,2) = d1notH(i+1,2) - d1notH(i,2);
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1aylH(i+1,1);
    data(i,2) = d1aylH(i+1,2) - d1aylH(i,2) + 1*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1cxyH(i+1,1);
    data(i,2) = d1cxyH(i+1,2) - d1cxyH(i,2)+2*k;
    data(i,3) = d1cxyH(i,3);
end
errorbar(data(:,1),data(:,2),data(:,3))
for i=18:-1:1
    data(i,1) = d1mgtH(i+1,1);
    data(i,2) = d1mgtH(i+1,2) - d1mgtH(i,2)+9*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1moqH(i+1,1);
    data(i,2) = d1moqH(i+1,2) - d1moqH(i,2)+3*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1qfmH(i+1,1);
    data(i,2) = d1qfmH(i+1,2) - d1qfmH(i,2)+4*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1qnfH(i+1,1);
    data(i,2) = d1qnfH(i+1,2) - d1qnfH(i,2)+5*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1qsaH(i+1,1);
    data(i,2) = d1qsaH(i+1,2) - d1qsaH(i,2)+6*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d1ygeH(i+1,1);
    data(i,2) = d1ygeH(i+1,2) - d1ygeH(i,2)+7*k;
end
plot(data(:,1),data(:,2))
for i=18:-1:1
    data(i,1) = d2driH(i+1,1);
    data(i,2) = d2driH(i+1,2) - d2driH(i,2)+8*k;
end
plot(data(:,1),data(:,2))
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