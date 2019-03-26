%Projekt PUST
%Zadanie 6
%Funkcja obliczaj?ca b??d DMC

clear all
lab1_zad3 = fullfile('Lab1Zad3b.m');

run(lab1_zad3);

%przypisanie odpowiedzi skokowej (znormalizowanej=
st = Ynorm(2:length(Ynorm));

addpath('F:\SerialCommunication'); % add a path to the functions
initSerialControl COM5 % initialise com port

% podstawowe wartosci
Upp = 28;
Ypp = 35;
iterNum = 700;
yZad = ones(iterNum, 1)*Ypp;
yZad(1:300) = 44;
yZad(301:700) = 39;
yZad = yZad - Ypp;
Umax=100;
Umin=0;

%REGULATOR DMC -----------------------------------------------------
%horyzonty
D = 734;
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

% Obliczanie parametrów regulatora
I=eye(Nu);
K=((M'*M+lambda*I)^(-1))*M';
Ku=K(1,:)*Mp;
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

ue = ke*e;
uu = Ku*dUpast;

du = ue-uu;
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
    
    sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
                     [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values

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