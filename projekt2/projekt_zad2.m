% obliczanie wzmocnienia statycznego torów wejście-wyjście i
% zakłócenie-wyjście

clear  all;

iterNum = 300;
licz_skok=10;

for i = 1 : licz_skok
    u = ones(iterNum, 1)*0*i;
    z = ones(iterNum, 1)*1*i;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
    figure(1);
    subplot(2,1,1);
    plot(y);
    xlabel('k');
    ylabel('y');
    ylim([0 12]);
    xlim([0 iterNum]);
    hold on;
    
    subplot(2,1,2);
    plot(z);
    xlabel('k');
    ylabel('z');
    xlim([0 iterNum]);
    ylim([-1 11]);
    hold on;
    
    yStat(i) = y(iterNum);
end
matlab2tikz('sprawko_dane/matlabtotikz/odp_skok_z.tex', 'showInfo', false);
hold off;



figure(2);
plot(1:licz_skok,yStat);
xlabel('Z statyczne');
ylabel('Y statyczne');
xlim([1 10]);
matlab2tikz('sprawko_dane/matlabtotikz/char_stat_Z_zmien.tex', 'showInfo', false);

wzmStatZ = (yStat(licz_skok)-yStat(1))/licz_skok;

for i = 1 : licz_skok
    u = ones(iterNum, 1)*1*i;
    z = ones(iterNum, 1)*0*i;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
    figure(3);
    subplot(2,1,1);
    plot(y);
    xlabel('k');
    ylabel('y');
    ylim([0 20]);
    xlim([0 iterNum]);
    hold on;
    
    subplot(2,1,2);
    plot(u);
    xlabel('k');
    ylabel('u');
    xlim([0 iterNum]);
    ylim([-1 11]);
    hold on;
    
    yStat(i) = y(iterNum);
end
matlab2tikz('sprawko_dane/matlabtotikz/odp_skok_u.tex', 'showInfo', false);

wzmStatU = (yStat(licz_skok)-yStat(1))/licz_skok;

hold off;
figure(4);
plot(yStat);
plot(1:licz_skok,yStat);
xlabel('U statyczne');
ylabel('Y statyczne');
xlim([1 10]);
matlab2tikz('sprawko_dane/matlabtotikz/char_stat_U_zmien.tex', 'showInfo', false);
