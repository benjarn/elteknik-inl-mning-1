%% Inlämning 1 - Elteknik 2014

%% A.a
%Fasströmmar
% 1.
U=400+j0;
phi1=acosd(0.98);
P1=250e3;

I1=P1/(sqrt(3)*U*cosd(phi1))

% 2.
S2=100e3;
phi2=acosd(0.95);

Ic=S2/(sqrt(3)*U);
I2=complement(Ic)

% 3.
phimotor=acosd(0.8);
Imotor=200;
Qbat=-90e3;

Pmotor=U*Imotor*cosd(phimotor);
Qmotor=U*Imotor*sind(phimotor);
Qtot=Qmotor+Qbat;
S3=Pmotor+jQtot;
I3=complement(S3/(sqrt(3)*U))

% Itot
Itot=I1+I2+I3

% Impedanser
Z=U/[I1 I2 I3 I4]

% Visardiagram
%figure()
% kanske med quiver?
% eller compass?

