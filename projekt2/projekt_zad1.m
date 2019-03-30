clear  all;

iterNum = 500;
skok=2;
u = zeros(iterNum, 1);
z = ones(iterNum, 1)*0;
y = ones(iterNum, 1)*skok;

for k = 7:iterNum
   y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
figure(1)
plot(y);
xlabel('k');
ylabel('y(k)');
title('sprawdzenie punktu pracy');

nazwa1 = sprintf('sprawko_dane/pkt_prac_z_u=%g.txt',skok);
 
file = fopen(nazwa1, 'w');
A = [(1:iterNum);y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);