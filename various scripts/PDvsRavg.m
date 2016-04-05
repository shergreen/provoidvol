x = csvread('packDensVSRprobe.txt');
%data that we want is standard deviation and average at each step size
count = zeros(1,16); avg = zeros(1,16);
moment2 = zeros(1,16);

for i=1:size(x,1)
    bin = 1+floor(5*x(i,2) + 0.5);
    avg(bin) = avg(bin) + x(i,1);
    moment2(bin) = moment2(bin) + x(i,1)*x(i,1);
    count(bin) = count(bin) + 1;
end

avg = avg ./ count;
moment2 = moment2 ./ count;
std = sqrt(moment2 - avg.*avg);

%plot(0.0:0.2:3.0,avg)
errorbar(0.0:0.2:3.0,avg,std)
xlabel('R_{probe}')
ylabel('Average Packing Density')