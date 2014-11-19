%% Inlamning 1 - Elteknik 2014

%% A.a
clear all;
close all;
clc;
cart2pol = @(z) [abs(z), angle(z)*180/pi]; % För att skriva ut värden på polär form
pol2cart = @(a, phi)(a*cosd(phi)+j*a*sind(phi));
%Fasstrommar
% 1.
U=400+j*0;
phi1=acosd(0.98);
P1=250e3;

I1r=P1/(sqrt(3)*U*cosd(phi1)); % realdel av I1
I1=pol2cart(I1r, -phi1)
cart2pol(I1)

% 2.
S2=100e3;
phi2=acosd(0.95);

I2r=S2/(sqrt(3)*U);
I2=pol2cart(I2r, -phi2)
cart2pol(I2)

% 3.
% Givna enheter
phimotor=acosd(0.8);
Imotor=200;
Qbat=-90e3;

Pmotor=sqrt(3)*U*Imotor*cosd(phimotor); % Räknar ut aktiv effekt för motorn
Qmotor=sqrt(3)*U*Imotor*sind(phimotor); % Reaktiv effekt för motorn
Qtot=Qmotor+Qbat; % Total reaktiv effekt
S3=Pmotor+j*Qtot; % Skenbar effekt
I3=conj(S3/(sqrt(3)*U))
cart2pol(I3)

% Itot
Itot=I1+I2+I3
cart2pol(Itot)

% Impedanser
Z=U./[I1 I2 I3 Itot]

% A.b
% Skenbar effekt för varje last samt den totala
Q1=U*real(I1)*sind(phi1);
S1=P1+j*Q1
S2=sqrt(3)*U*conj(I2)
S3
Stot=S1+S2+S3

%% B
clear all
close all
clc;
cart2pol = @(z) [abs(z), angle(z)*180/pi]

% Lägg in huvudspänningar med U1 som referens
U1=400+j*0;
U2=400*(cosd(-120)+j*sind(-120));
U3=400*(cosd(-240)+j*sind(-240));
% Impedanslaster
Z1=30+j*40;
Z2=100+j*0;
Z3=20-j*15;
% Ström genom verje impedans
Iab=U1/Z1
cart2pol(Iab)
Ibc=U2/Z2
cart2pol(Ibc)
Iac=U3/Z3
cart2pol(Iac)
% Fasströmmar
Ia=Iab-Iac
cart2pol(Ia)
Ib=Ibc-Iab
cart2pol(Ib)
Ic=Iac-Ibc
cart2pol(Ic)


%% C

%% D
clear all;
close all;
clc;
% D.a
% Givna storheter
S=1160;
U1=220;
U2=110;
% Märkströmmar
I1=S/U1
I2=S/U2
% D.b
bobbin=4; % mm
isolering=2;
distans=2;
% total bredd för b1+b2=30-(bobbin+isolering+distans)
b1b2=30-(bobbin+isolering+distans);
h=90-2*bobbin; % Effekttiv lindningshöjd
Al=b1b2*h; % Lindningsarea mm^2

% D.c
% intervall för Acu
Acu = Al.*[0.6, 0.7] % Koppararea

% D.d
Jp = [1.5, 1.8]; % Strömtäthet A/mm^2
Js = [1.9, 2.1];

A1=I1./Jp % ledningsarea för min och maxström
A2=I2./Js

d1=sqrt(A1./pi).*2 % Ledningsdiameter för min och maxström
d2=sqrt(A2./pi).*2

% D.e
% N1 = 2*N2
%N1*A1+N2*A2=Acu
N2(1)=floor(Acu(1)/(2*A1(1)+A2(1))); % N2 antal varv, avrundas neråt
N2(2)=floor(Acu(2)/(2*A1(2)+A2(2)));
N1=2.*N2
N2

% D.f
d1v = 2.00; % Vald diameter mm
d2v = 2.60;
N2v = sum(N2)/2; % Väljer antal varv i mitten av intervallet
N1v = 2*N2v;
AcuN1=N1v*pi*((d1v+0.1)/2)^2; % Koppararean för primärlindning
AcuN2=N2v*pi*((d1v+0.1)/2)^2; % -- sekundärlindning
Acuv=AcuN1+AcuN2; % totala koppararean
kcu=Acuv/Al % fyllfaktorn

b1v=AcuN1/(h*kcu) % b1 beräknad med fillfaktor och höjd
b2v=AcuN2/(h*kcu) % b2

% testar om de får plats
b1b2 - (b1v+b2v) >= 0
b1b2 - (b1v+b2v) % skriver ut värdet

% D.g
Bmax = [1.1, 1.2]; % Bmax min och max
f=50; % frekvens Hz
b=60; % bredden på järnkärnan inuti bobbin
kFe=[0.90, 0.95]; % kFe min och max
A=1e6*U1./(4.4*f*N1v*Bmax); % Area i mm
d=A/b; % effektiv järntjocklek
d/kFe(1) % total järntjocklek med fyllnadsfaktor 0.9
d/kFe(2) % -- 0.95

dv=77; % vald total tjocklek på bobbin
kFev=kFe(1);
Av=dv*b*kFev;
Bmaxv=1e6*U1./(4.4*f*N1v*Av)

% D.h
%Resistanser
rho=1.724e-5; % ohmmm^2/mm
r1=30+bobbin+b1v/2; % radie1 medelvärde
r2=30+bobbin+b2v/2;
l1=2*r1*pi*N1v; % längd på primärlindning, antar rund bobbin och medelvärdet på lindningarna
l2=2*r2*pi*N2v;
a1=pi*(d1v/2)^2; % ledarens koppararea
a2=pi*(d2v/2)^2;
R1=rho*l1/a1; % resistans på primärlindningen
R2=rho*l2/a2; % -- sekundärlindningen
Rk=R1+R2*2^2 % Rk=R1+R2'

% Läckreaktans
lm=(l1+l2)*1e-3/(N1v+N2v) ; % Medelvarvlängden m
delta=b1v/2+b2v/2+isolering; %avståndet mellan lindningarna mm
my0=4*pi*1e-7;
Xd=((2*pi*f*my0*lm*N1v^2)/h)*(delta+(b1v+b2v)/3)

% Kopparförluster
Pcu=Rk*(I2/2)^2

% Tomgångsförluster
Pfe=1.9; % W/kg
rho_fe=7.3; % kg/dm^3
Afe=(30+30+60+30+30)*(30+90+30)-2*(30*90); % Järnarea mm^2 (total - "fönster")
Vfe=Afe*dv*kFev; % järnvolym mm^3(Area*valt djup*fyllnadsfaktor)
Mfe=rho_fe*Vfe/1e6 % Järnmassa kg
P0=Pfe*Mfe % Tomgångsförlust

% Tomgångsström
%H*l=N2*I2;
H=170; % A/m
l=(2*30+2*90)*1e-3;
I02=(H*l)/N2v

% i
% Kortslutningstest
Uk=5.5;
Ik=5.3;
Pk=25;
Zk=Uk/Ik;
phik=acosd(Pk/(Ik*Uk));
Rkmeas=Zk*cosd(phik)
Xkmeas=Zk*sind(phik)


