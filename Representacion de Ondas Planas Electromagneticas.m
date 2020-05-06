% Simulación de Ondas Planas Electromagneticas
% Autor: José Rodrigo Fuentes Ramirez
% Versión: 1.0
clc; clear; close all;
% Caso er1 > er2
% Parámatros
n=100; %Número de muestras
f=1e9;
T=1/f;
c=3e8;
lambda = c/f;
thetaiengrados = 45;
thetai=degtorad(thetaiengrados);  
% Medio 1
er1=4;
mur1=1;
n1=sqrt(er1);
eta1=120*pi*sqrt(mur1/er1);
k1=(2*pi/lambda)*sqrt(mur1*er1);
% Medio 2
er2=1;
mur2=1;
n2=sqrt(er2);
eta2=120*pi*sqrt(mur2/er2);
k2 =(2*pi/lambda)*sqrt(mur2*er2);
%Hallamos el ángulo de reflexión y el ángulo de transmisión
thetar = thetai;
thetat = asin((n1/n2)*sin(thetai));
thetat_en_grados = radtodeg(thetat)
%Angulo Critico
thetac = asin(n2/n1);
theta_critico_grados = radtodeg(thetac)
%%%%%%%%ejes de coordenadas para los tres graficos%%%%%%%%%%%%%
Z1=-2*lambda:lambda/n:0; 
X1=0:lambda/n:2*lambda; 
[z1,x1]=meshgrid(Z1,X1);
X2 = lambda:lambda/n:3*lambda;
[z2,x2]=meshgrid(Z1,X2);
Z3 = 0:lambda/n:2*lambda;
[z3,x3]=meshgrid(Z3,X2);
%%%%%%%% Onda Incidente%%%%%%%%%%%%%
Ei0=1;
Eic =Ei0*exp(-1j*k1*(sin(thetai).*x1+cos(thetai).*z1));
%%figure(1)
subplot(2,2,3);
surf(z1,x1,real(Eic))
tituloOndaIncidente=sprintf('Onda Incidente,  Theta-inc = %dº ', thetaiengrados);
title(tituloOndaIncidente)
xlabel ('eje Z(m)'), ylabel ('eje X(m)'), zlabel('E(V)');
shading interp;
% Definimos el coeficiente de reflexión perpendicular
CoefRef = (eta2*cos(thetai)-eta1*cos(thetat))/(eta2*cos(thetai)+eta1*cos(thetat));
%Definimos la onda reflejada
Er0 = CoefRef*Ei0;
Erc = Er0*exp(-1j*k1*(sin(thetar).*x2-cos(thetar).*z2));
%Pintamos la onda reflejada
%%figure(2)
subplot(2,2,1);
surf(z2,x2,real(Erc))
tituloOndaReflejada=sprintf('Onda Reflejada,  Theta-ref = %dº ', thetaiengrados);
title(tituloOndaReflejada)
xlabel ('eje Z(m)'), ylabel ('eje X(m)'), zlabel('E(V)');
shading interp;
%Definimos el coeficiente de transmisión
CoefTra = (2*eta2*cos(thetai))/(eta2+cos(thetai)+eta1*cos(thetat));
%Definimos la onda transmitida
Et0 = CoefTra*Ei0;
Etc = Et0*exp(-1j*k2*(sin(thetat).*x3+cos(thetat).*z3)); 
%Campo total medio 1
 Etotal1 = Eic + Erc;
 %Pintamos el campo total
 %figure(5)
 subplot(2,2,4);
 surf(z1,x1,real(Etotal1))
 title('Campo total Medio 1')
 xlabel ('eje Z(m)'), ylabel ('eje X(m)'), zlabel('E(V)');
 shading interp;
% Si el angulo critico es menor que el incidente entonces no hay reflexion
% total
if(thetai < thetac)
    %Pintamos la onda Transmitida
    %figure(3)
    subplot(2,2,2);
    surf(z3,x3,real(Etc))
    tituloOndaTransmitida=sprintf('Onda Transmitida,  Theta-trans = %dº ', round(thetat_en_grados));
    title(tituloOndaTransmitida)
    xlabel ('eje Z(m)'), ylabel ('eje X(m)'), zlabel('E(V)');
    shading interp;
else
    %% Cuando tenemos reflexión total tenemos una onda de superficie en el medio 2
    alpha2 = k2*sqrt((er1/er2)*sin(thetai)^2 - 1);
    k2x = k2*sqrt(er1/er2)*sin(thetai);
    Etsup = Et0*exp(-alpha2.*z3-1j*k2x.*x3);   
    % Pintamos la onda superficial
    %figure(4)
    subplot(2,2,2);
    surf(z3,x3,real(Etsup))
    title('Onda Plana No Homogenéa en el medio 2')
    xlabel ('eje Z'), ylabel ('eje X'), zlabel('E(V)');
    shading interp;
 end

