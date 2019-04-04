%PUST Lab1
%Zad3b
%Wyznaczenie modelu odpowiedzi skokowej

%load('step_responseXX');
%load('step_response40.mat');


% 
% %ograniczenia nierownosciowe
% A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1];
% b = [100; 100; 100; 100; 0.1; 0.1; 0.1; 0.1];
% 
% %punkt startowy
% x0 = [1; 30; 5; 2];
% 
% %dodatkowe opcje
% options = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'iter');
% 
% %param = [T1 T2 K Td]
% x = fmincon(@model, x0, A, b, [], [], [], [], [], options);


%ograniczenia nierownosciowe
A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1];
b = [100; 100; 100; 100; -0.1; -0.1; -0.1; -1];

%punkt startowy
x0 = [76; 11; 0.1; 6];


options = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'off');

x = fmincon(@model, x0, A, b, [], [], [], [], [], options);

%dodatkowe opcje
%options = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 100);
%Td ca�kowite:
%IntCon = 4;

%param = [T1 T2 K Td]
%x = ga(@model, 4, A, b, [], [], [], [], [], IntCon, options);

x(4) = floor(x(4));


load('step-response80.mat');

    Ypp = 35.62;

    wartosc_skoku_U = 52; %40 - 28

    Ynorm = (step_response - Ypp)./wartosc_skoku_U;

    iterNum = length(step_response); %length(step_response)
    
    u = ones(1, iterNum);
    y = zeros(1,iterNum);

    alpha1 = exp(-1/x(1));
    alpha2 = exp (-1/x(2));
    a1 = - alpha1 - alpha2;
    a2 = alpha1 * alpha2;
    b1 = x(3) / (x(1)-x(2)) * (x(1) * (1-alpha1) - x(2)*(1-alpha2));
    b2 = x(3) / (x(1)-x(2)) * (alpha1*x(2) * (1-alpha2) - alpha2*x(1) * (1-alpha1));
    
    for k = 3+x(4):iterNum
        y(k) = b1 * u(k - x(4) - 1) + b2 * u(k - x(4) - 2) - a1 * y(k-1) - a2 * y(k-2);

    end
    
 figure(2)
 hold on
      stairs(y);
      hold off;
      hold on;
      plot(Ynorm);
      hold off;
 %Ynorm=y;
 
 
 nazwa1 = 'lab_odp_rzecz.txt';
 nazwa2 ='lab_odp_mod.txt';
 
 file = fopen(nazwa1, 'w');
 A = [(1:993);y];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);
 
file = fopen(nazwa2, 'w');
B = [(1:993);Ynorm];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);
      
function error = model(param)

    load('step-response80.mat');

    param(4) = floor(param(4));
    Ypp = 35.62;
    %wartosc_skoku_U = zalezy;
    wartosc_skoku_U = 52; %40 - 28

    Ynorm = (step_response - Ypp)./wartosc_skoku_U;

    iterNum = length(step_response); %length(step_response)
    %uzad = 1; %zadane sterowanie

    u = ones(1, iterNum);
    y = zeros(1,iterNum);

    alpha1 = exp(-1/param(1));
    alpha2 = exp (-1/param(2));
    a1 = - alpha1 - alpha2;
    a2 = alpha1 * alpha2;
    b1 = param(3) / (param(1)-param(2)) * (param(1) * (1-alpha1) - param(2)*(1-alpha2));
    b2 = param(3) / (param(1)-param(2)) * (alpha1*param(2) * (1-alpha2) - alpha2*param(1) * (1-alpha1));
    
    for k = 3+param(4):iterNum
        %y(k) = b1 * u(k - param(4) - 1) + b2 * u(k - param(4) - 2) - a1 * y(k-1) - a2 * y(k-2);
        y(k) = b1 * u(k - param(4) - 1) + b2 * u(k - param(4) - 2) - a1 * y(k-1) - a2 * y(k-2);

    end
    
    error = sum((Ynorm - y).^2);
    %disp('Wskaznik dopasowania: '+ error);
end

%dorobic optymalizacje i przedstawienie na wykresie