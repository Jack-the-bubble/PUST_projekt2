clear all

%pobranie odpowiedzi skokowych
run('projekt_zad3');
%wektor zawierający odpowiedź skokową zakłócenia
sz = yz(2:end);
%wektor zawierający odpowiedź skokową obiektu
s = y(2:end);
clear y u;

%czas symulacji
czas_sym=400;
%moment zadania wartości zadanej równej 1
moment_zad = 10;

%czy zakłócenie ma wystapić (0-nie, 1-tak)
zak=1;
%czy zakłócenie jest mierzone (0-nie, 1-tak)
pomiar_zak=1;
%moment wystąpienia zakłócenia
moment_zak = 134;
%czy zakłócenie ma mieć postać sinusoidy (0-nie, 1-tak)
zak_sin = 0;
%czy ma wystąpić dodatkowe zakłócenie w postaci szumu
szum = 1;
%współczynnik skalujacy wartości szumu
wsp_skal_szum = 0.5;
%amplituda zakłócenia w postaci sinusoidy
sin_A = 0.1;

%horyzont dynamiki regulatora
D=230;
%horyzont zakłócenia
Dz=229;%229;
%horyzont predykcji
N=50;
%horyzont sterowania
Nu=4;
%wspóczynnik kary za przyrosty sterowania
lambda=4;

%najlepsze - Marcin
% D=230;
% Dz=229;
% N=23;
% Nu=1;
% lambda=1;

%macierz M
M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end
   end
end

%macierz Mp
MP=zeros(N,D-1);
for i=1:N
   for j=1:D-1
      if i+j<=D
         MP(i,j)=s(i+j)-s(j);
      else
         MP(i,j)=s(D)-s(j);
      end      
   end
end

%macierz Mzp
MZP=zeros(N,Dz-1);
for i=1:N
   for j=1:Dz-1
      if i+j<=Dz
         MZP(i,j)=sz(i+j)-sz(j);
      else
         MZP(i,j)=sz(Dz)-sz(j);
      end      
   end
end
MZP = [sz(1:N) MZP];

%obliczanie parametrow regulatora: skalara ke, kz i wektora Ku
I=eye(Nu);
K=((M'*M+lambda*I)^(-1))*M';
ku=K(1,:)*MP;
kz=K(1,:)*MZP;
ke=sum(K(1,:));

%wektor przeszłych przyrostów sterowań
deltaup=zeros(D-1, 1);
%wektor przeszłych przyrostów wartości zakłócenia
deltazp=zeros(Dz, 1);
%sygnał wyjściowy obiektu
y = zeros(czas_sym, 1);
%trajektoria zadana
yzad=ones(czas_sym, 1)*1;
yzad(1:moment_zad)=0;
%sygnał sterujący
u = zeros(czas_sym, 1);
%aktualna zmiana zakłócenia
deltaz=0;
%aktualna zmiana sterowania
deltauk=0;

zaklocenie=zeros(czas_sym, 1);

%utworzenie odpowiedniego wektora z zakłóceniami
if zak_sin==0
    if zak==0
        zaklocenie=zeros(czas_sym, 1);
    elseif zak==1
        zaklocenie=ones(czas_sym, 1);
        zaklocenie(1:moment_zak)=0;
    end
elseif zak_sin == 1
    zaklocenie(1:moment_zak) = 0;
    for i=135:czas_sym
        zaklocenie(i)=sin_A*sin(i*pi/16);
    end
end
if szum == 1
   for i = 1: czas_sym
       zaklocenie(i)=zaklocenie(i)+randn()*wsp_skal_szum;
   end
end

%początek pętli
for i=7:czas_sym
   
   %pobranie wyjscia obiektu
   y(i) = symulacja_obiektu3y(u(i-5), u(i-6), zaklocenie(i-2), zaklocenie(i-3), y(i-1), y(i-2));
   
   %wyliczenie uchybu
   ek=yzad(i)-y(i);
   
   %wyliczenie wartosci sterowania przy różnych warunkach
   
   if pomiar_zak==1
      deltaz=zaklocenie(i)-zaklocenie(i-1);
      deltazp = [deltaz; deltazp(1:end-1)];
   end      
   
   deltauk=ke*ek-ku*deltaup;
   if pomiar_zak==1
      deltauk=deltauk-kz*deltazp;
   end         
   
   deltaup = [deltauk; deltaup(1:end-1)];
   u(i)=u(i-1)+deltauk; 
end

%wyliczenie wskaźnika jakości
wskaznikDMC = sum((yzad - y).^2);
disp(wskaznikDMC)

% Graficzna prezentacja wynik�w oblicze�
figure(4)   
subplot(2,1,1);
plot(y);
hold on;
plot(yzad);
hold off;
title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
legend('y','y_{zad}')
xlabel('k');
ylabel('y,y_{zad}');
subplot(2,1,2);
stairs(u);
xlabel('k');
legend('u');
ylabel('u');

if(zak==0)
    nazwa1 = sprintf('sprawko_dane/DMC_bez_zak/U__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,wskaznikDMC);
    nazwa2 = sprintf('sprawko_dane/DMC_bez_zak/Y__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,wskaznikDMC);
    nazwa3 = 'sprawko_dane/DMC_bez_zak/Yzad.txt';
else
    if(pomiar_zak == 1)
        if(zak_sin == 1)
            nazwa1 = sprintf('sprawko_dane/DMC_zak/U__DMC_E=%g_Dz=%g_sin=%g_.txt',wskaznikDMC,Dz,sin_A);
            nazwa2 = sprintf('sprawko_dane/DMC_zak/Y__DMC_E=%g_Dz=%g_sin=%g_.txt',wskaznikDMC,Dz,sin_A);
            nazwa3 = 'sprawko_dane/DMC_zak/Yzad.txt';
        else
            if(szum == 1)
                nazwa1 = sprintf('sprawko_dane/DMC_zak/U__DMC_E=%g_Dz=%g_szum=%g.txt',wskaznikDMC,Dz,wsp_skal_szum);
                nazwa2 = sprintf('sprawko_dane/DMC_zak/Y__DMC_E=%g_Dz=%g_szum=%g.txt',wskaznikDMC,Dz,wsp_skal_szum);
                nazwa3 = 'sprawko_dane/DMC_zak/Yzad.txt';    
            else    
                nazwa1 = sprintf('sprawko_dane/DMC_zak/U__DMC_E=%g_Dz=%g.txt',wskaznikDMC,Dz);
                nazwa2 = sprintf('sprawko_dane/DMC_zak/Y__DMC_E=%g_Dz=%g.txt',wskaznikDMC,Dz);
                nazwa3 = 'sprawko_dane/DMC_zak/Yzad.txt'; 
            end
        end
    elseif (pomiar_zak == 0)
        if(zak_sin == 0)
            nazwa1 = sprintf('sprawko_dane/DMC_zak/U__DMC_E=%g_bez_pom.txt',wskaznikDMC);
            nazwa2 = sprintf('sprawko_dane/DMC_zak/Y__DMC_E=%g_bez_pom.txt',wskaznikDMC);
            nazwa3 = 'sprawko_dane/DMC_zak/Yzad.txt';
        else
            nazwa1 = sprintf('sprawko_dane/DMC_zak/U__DMC_E=%g_bez_pom_sin=%g_.txt',wskaznikDMC,sin_A);
            nazwa2 = sprintf('sprawko_dane/DMC_zak/Y__DMC_E=%g_bez_pom_sin=%g_.txt',wskaznikDMC,sin_A);
            nazwa3 = 'sprawko_dane/DMC_zak/Yzad.txt';
        end
    end
end
 
% file = fopen(nazwa1, 'w');
% A = [(1:czas_sym);u'];
% fprintf(file, '%4.3f %.3f \n',A);
% fclose(file);
% 
% file = fopen(nazwa2, 'w');
% B = [(1:czas_sym);y'];
% fprintf(file, '%4.3f %.3f \n',B);
% fclose(file);
% 
% file = fopen(nazwa3, 'w');
% C = [(1:czas_sym);yzad'];
% fprintf(file, '%4.3f %.3f \n',C);
% fclose(file);
