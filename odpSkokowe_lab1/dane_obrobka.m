load('step-response28z20.mat');

file = fopen('odp_skok_dz=20.txt', 'w');
A = [(1:length(step_response));step_response];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen('z20.txt', 'w');
B = [(1:length(step_response));ones(length(step_response),1)'*20];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

clear all

load('step-response28z40.mat');

file = fopen('odp_skok_dz=40.txt', 'w');
A = [(1:length(step_response));step_response];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen('z40.txt', 'w');
B = [(1:length(step_response));ones(length(step_response),1)'*40];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

clear all

load('step-response28z60.mat');

file = fopen('odp_skok_dz=60.txt', 'w');
A = [(1:length(step_response));step_response];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen('z60.txt', 'w');
B = [(1:length(step_response));ones(length(step_response),1)'*60];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

clear all

file = fopen('char_stat_Y.txt', 'w');
A = [(1:4);[34 38 42.25 46]];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen('char_stat_U.txt', 'w');
B = [(1:4);[0 20 40 60]];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

