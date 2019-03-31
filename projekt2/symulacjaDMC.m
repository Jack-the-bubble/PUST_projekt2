print = 1;
to_save =1 ; % Czy zapisac


% Symulacja algorytmu DMC z mozliwoscia uwzgledniania zaklocen
% oraz szumow pomiarowych
Upp=0; Ypp=0; Zpp=0;
Tp = 0.5;

% Dostrajalne parametry regulatora
Dz = 230;          % Horyzont dynamiki zaklocen
D =230;          % Horyzont dynamiki obiektu
N = D;          % Horyzont predykcji
Nu = D;         % Horyzont sterowania
lambda = 1;       % Ograniczenie zmian sterowania

% Przygotowanie do testow
EOS = 300;
Yzad = zeros(EOS, 1);
Yzad(20:end) = 1;
Y = zeros(EOS, 1);
U = zeros(EOS, 1);

% Wszystko co zwiazane z zakloceniami
pomiar_zaklocenia =0;
skok_zaklocenia = 120; % Chwila aktywacji zaklocenia
Z = zeros(EOS, 1);
Z(skok_zaklocenia:end) = 0;  % Zwykle zaklocenie
%Zaklocenie sinusoidalne
% Tf=15;
% A=1;
% Z(skok_zaklocenia:end) =A*sin(linspace(0, 2*pi*(EOS - skok_zaklocenia)*Tp/Tf, EOS - skok_zaklocenia + 1)); 

% Szumy przy pomiarze zak��cenia
szumy_pomiarowe = 0;
Zz = Z;
Amp = 0.8;
jednorodny = 0;
% Zz = Zz + Amp*(rand(length(Z), 1) - 0.5); % Rozklad jednostajny
Zz = Zz + normrnd(0, Amp/4, length(Z), 1);  % Rozklad normalny

% Wczytywanie odpowiedzi skokowej
S = y(2:end);
Sz = yz(2:end);

% Tworzenie macierzy M, MP i ZP
M = zeros(N, Nu);
for i = 1:N
    for j = 1:Nu
        if(j<=i)
            M(i,j) = S(i-j+1);
        end
    end
end

MP = zeros(N, D-1);
for i = 1:D-1
    for j = 1:N
        if (i + j) <= D-1
            MP(j,i) = S(i+j) - S(i);
        else
            MP(j,i) = S(D) - S(i);
        end
    end
end

MzP = zeros(N, Dz);
for i = 1:Dz
    for j = 1:N
        if i == 1
            if j <= Dz
                MzP(j,i) = Sz(i);
            else
                MzP(j,i) = Sz(Dz);
            end
        elseif (i-1 + j) <= Dz-1
            MzP(j,i) = Sz(i-1+j) - Sz(i-1);
        else
            MzP(j,i) = Sz(Dz) - Sz(i-1);
        end
    end
end

% Wektor przeszlych przyrostow sterowania
deltaUp = zeros(D-1, 1);   
deltaZp = zeros(Dz, 1);

% Wartosci ktore moza obliczyc przed symulacja:
K = ((M'*M + lambda*eye(Nu))^-1)*M';
ke = sum(K(1,:));
ku = K(1,:)*MP;
kz = K(1,:)*MzP;
for k = 7:EOS
    % Wyjscie obiektu
    Y(k)=symulacja_obiektu3y(U(k-5), U(k-6), Z(k-2), Z(k-3), Y(k-1), Y(k-2));
    
    % Zaklocenie
    if szumy_pomiarowe == 0
        deltaZk = Z(k) - Z(k-1);
    else
        deltaZk = Zz(k) - Zz(k-1);
    end
    deltaZp = [deltaZk; deltaZp(1:end-1)];
    
    % Algorytm DMC
    if pomiar_zaklocenia == 0
        deltaUk = ke*(Yzad(k) - Y(k)) - ku*deltaUp; 
    else
        deltaUk = ke*(Yzad(k) - Y(k)) - ku*deltaUp - kz*deltaZp; 
    end
    
    % Prawo regulacji
    U(k) = U(k-1) + deltaUk;

    % Zapis przeszlych przyrostow sterowania
    deltaUp = [deltaUk; deltaUp(1:end-1)];
end
 
format long g
E = sum((Yzad - Y).^2)



% E = sum(abs(Yzad - Y));
 
if print
    figure(1);
    subplot(2,1,1)
    plot(1:k, Y(1:k), 'LineWidth', 1.1); % Wyjscie obiektu
    hold on
    plot(1:k, Yzad(1:k), 'LineWidth', 1.1); % zadana
    hold off
    xlim([1, EOS]);
    title('Sygnal wyjsciowy');
    xlabel('Numer probki (k)');
    grid on;
    
    subplot(2,1,2)
    plot(1:k, U(1:k),'LineWidth', 1.1); % sterowanie
    hold on
    if szumy_pomiarowe == 0
        plot(1:k, Z(1:k),'LineWidth', 1.1); % zaklocenie
    else
        plot(1:k, Zz(1:k),'LineWidth', 1.1); % zaklocenie
    end
    hold off
    xlim([1, EOS]);
    title('Sygnal sterujacy oraz zaklocenie');
    xlabel('Numer probki (k)');
    grid on;
end
    
% 
% 
% 
% if (to_save == 1)
%     % Zapis do pliku
%     fileU = sprintf('dane/DMCszumy/DMCUszumy%.2fE%.4f.txt', Amp, E);
%     fileZ = sprintf('dane/DMCszumy/DMCZszumy%.2fE%.4f.txt', Amp, E);
%     fileY = sprintf('dane/DMCszumy/DMCYszumy%.2fE%.4f.txt', Amp, E);
%     fileYzad = sprintf('dane/DMCszumy/DMCYzadszumy%.2fE%.4f.txt', Amp, E);
%     fileIDU = fopen(fileU,'w');
%     fileIDZ = fopen(fileZ,'w');
%     fileIDY = fopen(fileY,'w');
%     fileIDYzad = fopen(fileYzad,'w');
% 
%     for j=1:length(U)
%         fprintf( fileIDU,    '%d  %f\r\n', j, U(j));
%         if szumy_pomiarowe == 0
%             fprintf( fileIDZ,    '%d  %f\r\n', j, Z(j));
%         else
%             fprintf( fileIDZ,    '%d  %f\r\n', j, Zz(j));
%         end
%         fprintf( fileIDY,    '%d  %f\r\n', j, Y(j));
%         fprintf( fileIDYzad, '%d  %f\r\n', j, Yzad(j));
%     end
%     fclose(fileIDU);
%     fclose(fileIDZ);
%     fclose(fileIDY);
%     fclose(fileIDYzad);
%     
%     % Gotowiec do pliku tex
%     fprintf('\n\\tikzsetnextfilename{}\n');
%     fprintf('\\begin{figure}[tb]\n');
%     fprintf('	\\centering\n')
%     fprintf('	\\tikzsetnextfilename{DMC/%s}\n', name);
%     fprintf('	\\begin{tikzpicture}\n');
%     fprintf('	\\begin{groupplot}[group style={group size=1 by 2,vertical sep=\\odstepionowy},\n');
%     fprintf('	width=\\szer,height=\\wys]\n');
%             %%1
%     fprintf('		\\nextgroupplot\n');
%     fprintf('		[xmin=0, xmax=%d, \n', EOS);
%     fprintf('		xlabel=$k$,ylabel={$y^{\\mathrm{zad}}, \\ y$},legend cell align=left, legend pos=south east,]\n');
% 
%     fprintf('		\\addplot[color=blue,semithick] file {../../../Projekt/Matlab/%s};\n',fileY);
%     fprintf('		\\addplot[const plot, color=red,semithick,densely dashed] file {../../../Projekt/Matlab/%s};\n', fileYzad);		
%     fprintf('		\\legend{$y$,$y^{\\mathrm{zad}}$}\n');
%             %%2
%     fprintf('		\\nextgroupplot\n');
%     fprintf('		[xmin=0,xmax=%d,xlabel=$k$,ylabel={$u, \\ z$},legend cell align=left, legend pos=south east,]\n', EOS);
%     fprintf('		\\addplot[const plot, color=black,semithick]\n');
%     fprintf('		file {../../../Projekt/Matlab/%s};\n', fileU);
%     fprintf('		\\addplot[color=blue,semithick]\n');
%     fprintf('		file {../../../Projekt/Matlab/%s};\n', fileZ);
%     fprintf('		\\legend{$u$,$z$}\n');
%     fprintf('	\\end{groupplot}\n');
%     fprintf('	\\end{tikzpicture}\n');
%     fprintf('\\end{figure}\n\n');
%     
%     % I wstawianie w spawku
%     fprintf('\n\\begin{figure}[tb]\n');
% 	fprintf('\\centering\n');
% 	fprintf('\\includegraphics[scale=1]{generacja_wykresow/zapisz_pdf/DMC/%s}\n', name);
%     if (jednorodny==1)
% 	fprintf('\\caption{Przebieg sygna�u wej�ciowego i wyj�ciowego przy dodaniu szum�w pomiarowych, dla sta�ego skoku zak��cenia o rozk�adzie $\\mathcal{U} \\sim (%.2f, \\ %.2f)$.\n ',-Amp/2, Amp/2);
%     else
% 	fprintf('\\caption{Przebieg sygna�u wej�ciowego i wyj�ciowego przy dodaniu szum�w pomiarowych, dla sta�ego skoku zak��cenia o rozk�adzie $\\mathcal{N} \\sim (0, \\ %.2f)$.\n ',Amp/4);       
%     end
%     fprintf('\\ Wska�nik jako�ci $E=\\num{%.4f}$.}\n', E);
% 	fprintf('\\label{%s}\n', name);
%     fprintf('\\end{figure}\n');
% end


