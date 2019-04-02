load('DMC_N=170_Nu=40_lambda=2_zak_niemierz_15.mat');
nazwa1 = 'Y__dmc_zak_niemierz_skok_15.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'U__dmc_zak_niemierz_skok_15.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'Yzad__dmc_zak_niemierz_skok_15.txt';
file = fopen(nazwa3, 'w');
A = [(1:length(Y));yZad'+Ypp];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

%----------------------------------------------------
load('DMC_N=170_Nu=40_lambda=2_zak_niemierz_30.mat');
nazwa1 = 'Y__dmc_zak_niemierz_skok_30.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'U__dmc_zak_niemierz_skok_30.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'Yzad__dmc_zak_niemierz_skok_30.txt';
file = fopen(nazwa3, 'w');
A = [(1:length(Y));yZad'+Ypp];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

%----------------------------------------------------
load('DMC_N=170_Nu=40_lambda=2_zak_mierz_15.mat');
nazwa1 = 'Y__dmc_zak_mierz_skok_15.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'U__dmc_zak_mierz_skok_15.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'Yzad__dmc_zak_mierz_skok_15.txt';
file = fopen(nazwa3, 'w');
A = [(1:length(Y));yZad'+Ypp];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

%----------------------------------------------------
load('DMC_N=170_Nu=40_lambda=2_zak_mierz_30.mat');
nazwa1 = 'Y__dmc_zak_mierz_skok_30.txt';
file = fopen(nazwa1, 'w');
A = [(1:length(Y));Y'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa2 = 'U__dmc_zak_mierz_skok_30.txt';
file = fopen(nazwa2, 'w');
A = [(1:length(Y));U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

nazwa3 = 'Yzad__dmc_zak_mierz_skok_30.txt';
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