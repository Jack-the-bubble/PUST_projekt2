% odpowiedź skokowa dla algorytmu DMC

clear  all;

iterNum = 300;

% skok zakłóceń

    u = ones(iterNum, 1)*0;
    z = ones(iterNum, 1)*1;
    yz = ones(iterNum, 1)*0;
for k = 7:iterNum
    yz(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), yz(k-1), yz(k-2));
end
%     figure(1);
%     
%     plot(yz)

    % skok sterowania

    u = ones(iterNum, 1)*1;
    z = ones(iterNum, 1)*0;
    y = ones(iterNum, 1)*0;
for k = 7:iterNum
    y(k) = symulacja_obiektu3y(u(k-5), u(k-6), z(k-2), z(k-3), y(k-1), y(k-2));
end
%     figure(2);
  %   plot(y)
    
% nazwa1 = 'sprawko_dane/odp_dmc_u.txt';
% nazwa2 = 'sprawko_dane/odp_dmc_z.txt';% 
% file = fopen(nazwa1, 'w');
% A = [(1:iterNum);y'];
% fprintf(file, '%4.3f %.3f \n',A);
% fclose(file);
% 
% file = fopen(nazwa2, 'w');
% B = [(1:iterNum);yz'];
% fprintf(file, '%4.3f %.3f \n',B);
% fclose(file);