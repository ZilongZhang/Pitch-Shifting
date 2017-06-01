clear id
clear res
close all
figure
hold on 
lx = 1024;
n4 = 256
for i = 1:10
id(i,:) = zeros(n4*8+lx,1);
id(i,n4*(i-1)+1:n4*(i-1)+lx) = hanning(lx)*0.5;
plot(id(i,:));
end
res(1,:) = id(1,:);
figure
plot(res(1,:));
hold on
for i = 2:10
    res(i,:) = id(i,:) + res(i-1,:);
    plot(res(i,:));
end

