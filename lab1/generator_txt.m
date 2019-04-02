load('DMC_N=150_Nu=15_lambda=01.mat');
nazwa1 = 'txt/Y__dmc_N=150_Nu=15_lambda=01.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'txt/U__dmc__N=150_Nu=15_lambda=01.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'txt/Yzad__dmc__N=150_Nu=15_lambda=01.txt';
file = fopen(nazwa3, 'w');
A = [(1:length(Y));yZad'+Ypp];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

%----------------------------------------------------
load('DMC_N=170_Nu=40_lambda=2.mat');
nazwa1 = 'txt/Y__dmc_N=170_Nu=40_lambda=2.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'txt/U__dmc_N=170_Nu=40_lambda=2.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'txt/Yzad__dmc_N=170_Nu=40_lambda=2.txt';
file = fopen(nazwa3, 'w');
A = [(1:length(Y));yZad'+Ypp];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

figure(17)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U);