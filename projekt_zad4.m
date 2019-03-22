clear all

run('projekt_zad3');

sz = yz;
s = y(2:end);

clear y u;


zad6 = 0;
RAND = 1;
% Horyzonty

% najlepsze bez zakłóceń
D=210;
Dz=216;
N=23;
Nu=1;



% Wsp�czynnik kary za przyrosty sterowania

lambda=1;

% Generacja macierzy

M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end;
   end;
end;

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

MZP=zeros(N,Dz-1);
for i=1:N
   for j=1:D-1
      if i+j<=D
         MZP(i,j)=sz(i+j)-sz(j);
      else
         MZP(i,j)=sz(D)-sz(j);
      end      
   end
end

% Obliczanie parametr�w regulatora

I=eye(Nu);
K=((M'*M+lambda*I)^-1)*M';
ku=K(1,:)*MP;
kz=K(1,:)*MZP;
ke=sum(K(1,:));


% % wyznaczanie warunków początkowych

czas_sym=400;%10

% Zak��cenie mierzalne
fl_pomiar_zak=1; % 1 - pomiar zak��cenia wykorzystany w regulatorze
% chwila_usterk=30;
% Warunki pocz�tkowe

yzad=ones(czas_sym, 1)*1;

if zad6==0
zaklocenie=ones(czas_sym, 1);
zaklocenie(1:134)=0;
elseif zad6 == 1
    zaklocenie(1:134) = 0;
    for i=135:czas_sym
        zaklocenie(i)=sin(i*pi/16);
    end
end

if RAND == 1
   for i = 1: czas_sym
       zaklocenie(i)=zaklocenie(i)+randn()*0.1;
   end
end

deltaup=zeros(D-1, 1);
deltazp=zeros(Dz-1, 1);

% G��wna p�tla programu
wyy = zeros(czas_sym, 1);
wyu = zeros(czas_sym, 1);
uk = 0;
for i=7:czas_sym
   
   % Obiekt
   wyy(i) = symulacja_obiektu3y(wyu(i-5), wyu(i-6), zaklocenie(i-2), zaklocenie(i-3), wyy(i-1), wyy(i-2));
   % Regulator DMC
   
   ek=yzad(i)-wyy(i);
   
   if fl_pomiar_zak==1
      for n=D-1:-1:2;
         deltazp(n)=deltazp(n-1);
      end
      deltazp(1)=zaklocenie(i)-zaklocenie(i-1);
   end      
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Prawo regulacji
   
   deltauk=ke*ek-ku*deltaup;
   if fl_pomiar_zak==1
      deltauk=deltauk-kz*deltazp;
   end         
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   deltaup = [deltauk; deltaup(1:end-1)];
   wyu(i)=wyu(i-1)+deltauk;
   
%    wyu(i)=uk;
end

wskaznikDMC = sum((yzad - wyy).^2);


% Graficzna prezentacja wynik�w oblicze�

figure(4)   
subplot(2,1,1);
plot(wyy);
hold on;
plot(yzad);
hold off;
title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
legend('y','y_{zad}')
xlabel('k');
ylabel('y,y_{zad}');
subplot(2,1,2);
stairs(wyu);
xlabel('k');
legend('u');
ylabel('u');


