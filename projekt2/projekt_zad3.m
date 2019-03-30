% odpowiedź skokowa dla algorytmu DMC

clear  all;

iterNum = 500;

% skok zakłóceń

    u = ones(iterNum, 1)*0;
    z = ones(iterNum, 1)*1;
    yz = ones(iterNum, 1)*0;
for k = 7:iterNum
    yz(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), yz(k-1), yz(k-2));
end
    figure(1);
    
    plot(yz)

    % skok sterowania

    u = ones(iterNum, 1)*1;
    z = ones(iterNum, 1)*0;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
    figure(2);
    plot(y)
    
    % nazwa1 = sprintf('dane_zad_5/DMC/U__DMC_D=%g_N=%g_Nu=%g_L=%g.txt',D,N,Nu,lambda);
% nazwa2 = sprintf('dane_zad_5/DMC/Y__DMC_D=%g_N=%g_Nu=%g_L=%g.txt',D,N,Nu,lambda);
% nazwa3 = 'dane_zad_5/DMC/Yzad.txt';
% 
% file = fopen(nazwa1, 'w');
% A = [(1:iterNum);U'];
% fprintf(file, '%4.3f %.3f \n',A);
% fclose(file);
% 
% file = fopen(nazwa2, 'w');
% B = [(1:iterNum);Y'];
% fprintf(file, '%4.3f %.3f \n',B);
% fclose(file);
% 
% file = fopen(nazwa3, 'w');
% C = [(1:iterNum);(yZad+Ypp)'];
% fprintf(file, '%4.3f %.3f \n',C);
% fclose(file);