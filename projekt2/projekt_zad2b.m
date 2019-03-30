clear  all;

iterNum = 300;
licz_skok=10;
% wykres charakterystyki statycznej Y(u, z)

for i = 1 : licz_skok
    u = ones(iterNum, 1)*1*i;
    
    y = ones(iterNum, 1)*0;
    for j = 1 : licz_skok
        z = ones(iterNum, 1)*1*j;
        for k = 7:iterNum
            y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
        end
        yStat(i, j) = y(iterNum);
    end
    
end
hold off;
figure(1);
surf(yStat);
xlabel('z');
ylabel('u');
zlabel('y(u, z)');
matlab2tikz('sprawko_dane/matlabtotikz/char_stat_X_Z_zmien.tex', 'showInfo', false);
