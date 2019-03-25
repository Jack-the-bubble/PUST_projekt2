% obliczanie wzmocnienia statycznego torów wejście-wyjście i
% zakłócenie-wyjście

clear  all;

iterNum = 500;


for i = 1 : 40
    u = ones(iterNum, 1)*0*i;
    z = ones(iterNum, 1)*1*i;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
    figure(1);
    plot(y)
    hold on;
    yStat(i) = y(iterNum);
end
hold off;
figure(2);
plot(yStat);

wzmStatZ = (yStat(40)-yStat(1))/40;

for i = 1 : 40
    u = ones(iterNum, 1)*1*i;
    z = ones(iterNum, 1)*0*i;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
    figure(3);
    plot(y)
    hold on;
    yStat(i) = y(iterNum);
end

wzmStatU = (yStat(40)-yStat(1))/40;

hold off;
figure(4);
plot(yStat);