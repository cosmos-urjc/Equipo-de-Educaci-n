%% Parte 1: Matlab como calculadora 
% Algo de teor�a 
% Las operaciones se ejecutan exactamente como en cualquier calculadora cient�fica, declarando unos datos (en unos "contenedores" llamados variables) 
% y operando sobre ellos. Algunos ejemplos: 

% Scritps vs Command Window: FIGHT!
% Breakpoints y secciones. Use them!
% Variables: nuestros contenedores. Nombres y palabras reservadas.
% Errores. Aprovechemos el int�rprete!
% Un espacio de trabajo limpio: clc y clear.

% Esto es un primer comentario. 

% Aritm�tica modular 
variable_1 = 1;                     %Primera variable de prueba �Qu� pasa si no ponemos ;? %�Por qu� comentar el c�digo y como hacerlo?
variable_2 = 2;                     %Segunda variable de prueba 
%Declara aqu� una tercera variable
variable_3 = 10e3;

sum = variable_1+variable_2;       %Suma y guarda el resultado en sum
dif = variable_1-variable_2;       %Resta y guarda el resultado en diff
pro = variable_1*variable_2;       %Multiplica y guarda el resultado en pro 
div = variable_1/variable_2;       %Divide y guarda en div

%Con jerarqu�a de operaciones (utilizando par�ntesis), calcula primero la suma de sum y dif con la tercera variable, 
%eleva al cuadrado el resultado y despu�s hazle la ra�z quinta.
result = ((sum+dif+variable_3)^2)^(1/5);
disp('Vamos a calcular el resultado!');
fprintf("El resultado de las anteriores operaciones es %.f \n", result);
fprintf("El resultado ha sido correcto");

%%
%Operaciones m�s avanzadas
pot = variable_1^variable_2;       %Potencias 
sq = sqrt(variable_1);             %Ra�z cuadrada (solo la cuadrada!) 
variable_1 = -1;                   %Reasignaci�n de valor a la variable variable_1
sqC = sqrt(variable_1);            %Ra�z cuadrada (en el plano complejo!)

%%
%�C�mo se calculan las funciones seno o coseno? ----> https://es.mathworks.com/

% Operaciones con matrices (the best!)
% Para MATLAB, todos los objetos que se manejan (pr�cticamente) son  tratados como arrays multidimensionales 
% (un escalar es una matriz de 1x1, un vector un array de nx1 o 1xn...). Esto hace que las operaciones matriciales 
% sean mucho m�s sencillas y r�pidas.

%Definici�n de arrays
scalar = 2;
vector = [1; 1; 2];
matrix = [1 2 3 4; 8 0 0 0];

matrix_1 = rand(scalar, scalar);
matrix_2 = eye(scalar);
matrix_3 = magic(scalar);

%Indexaci�n de arrays 
prueba = vector(1); 
prueba_2 = vector(1:2);
prueba_3 = matrix_1(:,1); 
prueba_4 = matrix_2(1,1);

concatenacion = [matrix_1 matrix_2; matrix_3 zeros(scalar, scalar)];

%Opera ahora
%ExProd(1:3,1) = scalar*vector;                 %Calcula el producto de un vector y un escalar
%ExProd(:,2) = ExProd(1,:)-vector;              %Algunas operaciones

%InProd(1) = dot(ExProd(1:2,2), vector(1:2,1));      %Producto escalar entre dos vectores
%InProd(2) = dot(ExProd(2,:), ExProd(2,:));          %Producto escalar del mismo vector consigo mismo (distancia eucl�dea!)
%InProd(3) = norm(ExProd(2,:))^2;                    %Producto escalar del mismo vector consigo mismo (distancia eucl�dea!)

MatrixProd = matrix_1*matrix_2*matrix_3;            %Producto de matrices
Error = matrix_1*vector(1:2,1);                     %Comprueba las dimensiones!

%%
%Cosas chulas de �lgebra 
vector = vector.';                                  %Traspuesta de un vector/matriz
vector = vector.^scalar;                            %Operaci�n vectorizada (componente a componente)
Det = det(matrix_1);                                %Determinante de una matriz
InvMatrix1 = inv(matrix_1);                         %Inversa de una matriz 
InvMatrix2 = matrix_1^(-1);                         %Inversa de una matriz. Y si ponemos .^? Observa la diferencia en el c�mputo num�rico!
[W, eigs] = eig(InvMatrix2);                        %An�lisis del espectro de una matrix

%%
%Definiciones importantes 
xmax = 1;                                           %Longitud m�xima
dx = 1e-3;                                          %Paso 
vector_1 = linspace(-5, 5, 100);                      %Divisi�n regular en N pasos del intervalo -5 a 5
vector_2 = 0:dx:xmax;                               %Divisi�n regular
vector_3 = [1 2 3];                                 %Vector

[X1, X2, X3] = meshgrid(vector_1, vector_2, vector_3);  %Producto cartesiano o mallado del espacio

% Algo de polinomios
pol = [1 2 3]; 
pol2 = poly(pol);
solutions = roots(pol2);

%%
% Ploteo y gr�ficas: https://www.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html
x = 0:1e-1:1;
y1 = 2+4*x; 
y2 = 2+3*x+6*x.^2;

figure(1) 
hold on 
subplot(2,1,1)
plot(x,y1, 'b'); 
subplot(2,1,2)
plot(x,y2, '.-r');
hold off 
legend('Recta', 'Par�bola'); 
xlabel('Eje x');
ylabel('Eje y');
title('Gr�fica de prueba'); 
grid on;

%%
% Ejercicio 0
% Sergio tiene 3 gatos y dos perros. En total, pesan 30 kg. Los dos perros menos un gato
% pesan 10 kg. Todos los perros pesan lo mismo, igual para los gatos. Resuleve para ambos pesos en forma matricial.
A = [3 2; -1 2];
b = [30; 10];

% Si ahora dos gatos y cuatro perros pesan 7.5 kilos, calcula la soluci�n de m�nimos cuadrados.
A = [A; 2 4]; 
b = [b; 7.5]; 
sol = pinv(A)*b;

%%
% Ejercicio 1
% Representar la trayectoria de un avi�n que describe una h�lice de radio 10, ascendiendo por el eje z.
% PISTA: para pasar 3 coordenadas, utiliza la funci�n plot3 en vez de plot.
t = linspace(0,2*pi,100);
z = linspace(0,10,100);
x = sin(t);               
y = cos(t);
figure(1)
view(3)
hold on
plot3(x,y,z); 
plot3(x(50), y(50), z(20), 'or');
hold off
grid on; 

%%
% Ejercicio 2
% Vamos a calcular el campo de velocidades de un determinado flujo
% potencial (aerodin�mico, el�ctrico, gravitatorio (?)). El espacio de configuraci�n es de 2 x 2 metros (R2). 
% Conocemos la posici�n de tres fuentes (una carga positiva), de dos sumideros (carga negativa), en (cos(45�), -sin(45�))
% y un v�rtice. Las cargas positivas son de 1,2 y 3 u; las negativas, -1, -3, -5; el v�rtice posee vorticidad 10.

% Posici�n de las cargas: (0,1), (1,0), (0,-1)
% Posici�n de los sumideros: (cos(45�), -sin(45�)), (-1.5, -1)
% Posici�n de un v�rtice: (0,0)

% Campo generado por una carga positiva: (qi/4*pi) * (1/||r-r_0||^2)
% Campo generado por una carga negativa: (qi/4*pi) * (1/||r-r_0||^2)
% Campo generado por un v�rtice:    (gamma*i/2*pi) * (log z-z_0)

% Obt�n la representaci�n gr�fica de las l�neas de campo. 
L = 2;                      %L�mite del espacio de configuracion
h = 1e-3;                   %Paso espacial
x = -L:h:L-h;               %Dominio horizontal
y = -L:h:L-h;               %Dominio vertical
[X, Y] = meshgrid(x,y);     %Mallado
Z = X-1i*Y;                 %Isomorfismo entre R2 y C

qp = [1 2 3]; 
qn = -[3 5]; 
Gamma = 100;

Rp = [1*1i; 1; -1*1i];
Rn = [cos(deg2rad(45))-1i*sin(deg2rad(45)); -1.5-1*1i];
Rg = 0;

F_v = ((1i*Gamma)/(2*pi))*(log(Z-Rg));                                                                      %Potencial torbellino
F_fp = ((qp(1))/(2*pi))*(log(Z-Rp(1)))+((qp(2))/(2*pi))*(log(Z-Rp(2)))+((qp(3))/(2*pi))*(log(Z-Rp(3)));     %Potencial fuente
F_fn = ((qn(1))/(2*pi))*(log(Z-Rn(1)))+((qn(2))/(2*pi))*(log(Z-Rn(2)));                                     %Potencial fuente
F_t = F_v+F_fp+F_fn;

phi = real(F_t);                           %Definici�n de funci�n potencial
psi = imag(F_t);                           %Definici�n de funci�n de corriente

%Representaci�n
figure(1)
plot3(X,Y,phi);  

figure(2) 
dens = 20;
streamslice(X,Y,phi,psi,dens);

figure(3) 
hold on
contour(X,Y,phi,dens);
contour(X,Y,psi,dens);
hold off