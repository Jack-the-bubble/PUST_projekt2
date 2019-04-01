%Projekt PUST
%Zadanie 6
%Funkcja obliczaj?ca b??d DMC

clear all
lab1_zad3 = fullfile('Lab1Zad3b.m');
lab2_zakl = fullfile('aproksymacja_odp_skok.m');

run(lab2_zakl);
run(lab1_zad3);

%przypisanie odpowiedzi skokowej (znormalizowanej=
st = Ynorm(2:length(Ynorm));
sz = YnormZ(2:length(YnormZ));

addpath('F:\SerialCommunication'); % add a path to the functions
initSerialControl COM5 % initialise com port

% podstawowe wartosci
Upp = 28;
Ypp = 34;
iterNum = 700;
yZad = ones(iterNum, 1)*Ypp;
yZad(1:700) = 44;
%yZad(301:700) = 39;
yZad = yZad - Ypp;
Umax=100;
Umin=0;
zak = 40;
pomiar_zak = 0;
a = 30; %amplituda zaklocenia

%REGULATOR DMC -----------------------------------------------------
%horyzonty
D = 734;
Dz = 
N = 170;
Nu = 40;
lambda = 2;

%PARAMETRY 
du = 0;
upast = 0.0; %poprzednia wartosc sterowania
e = 0.0; %uchyb

 u = zeros(iterNum, 1);
 U = ones(iterNum, 1)*Upp; %zmienic tak, zeby bral od poczatku aktualne wartosci u i y
 Y = ones(iterNum, 1)*Ypp;
dUpast = zeros(D-1, 1); %wektor przeszlych przyrostow sterowan
%wektor przeszłych przyrostów wartości zakłócenia
dZpast = zeros(Dz, 1);
%aktualna zmiana zakłócenia
dz=0;

zaklocenie=zeros(czas_sym, 1);
zaklocenie(250:end) = a;

% Macierz M
M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=st(i-j+1);
      end
   end
end

% Macierz Mp
Mp=zeros(N,D-1);
for i=1:N
   for j=1:D-1
      if (i+j)<=D-1
         Mp(i,j)=st(i+j)-st(j);
      else
         Mp(i,j)=st(D)-st(j);
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

% Obliczanie parametr�w regulatora
I=eye(Nu);
K=((M'*M+lambda*I)^(-1))*M';
Ku=K(1,:)*Mp;
kz=K(1,:)*MZP;
ke=sum(K(1,:));

% -------------- DO REGULACJI ---------------
k=2;
while(1)
upast = u(k-1);

%pobranie wyjscia obiektu
measurements = readMeasurements(1:7);
Y(k)= measurements(1);
y = Y(k)-Ypp;

e = yZad(k) - y;

if pomiar_zak==1
      dz=zaklocenie(i)-zaklocenie(i-1);
      dZpast = [dz; dZpast(1:end-1)];
   end

ue = ke*e;
uu = Ku*dUpast;

du = ue-uu;

if pomiar_zak==1
      du=du-kz*dZpast;
end

u(k) = upast+du;
U(k) = u(k)+Upp;

if U(k) <  Umin 
     U(k) = Umin;
     du = Umin-U(k-1);

elseif U(k) > Umax 
     U(k) = Umax;
     du = Umax - U(k-1);

end

    dUpast = [du; dUpast(1:end-1)];
    
    
    if k < 250
    sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
                     [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values
    else
    sendControlsToG1AndDisturbance(U(k),zaklocenie(k));
    
    end
    
%wykresy
        figure(17)
        subplot(2,1,1);
        plot(Y(1:k));
        hold on;
        plot(yZad(1:k)+Ypp);
        hold off;
        title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U(1:k));
        drawnow;
        
        disp("U: " + U(k) +" Y: "+ Y(k) +" Yzad: "+ yZad(k)+Ypp);
    
        k = k+1;
        waitForNewIteration();


end