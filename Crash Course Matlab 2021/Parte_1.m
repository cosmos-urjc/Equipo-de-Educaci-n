%% Parte 1: Matlab como calculadora 
% Algo de teoría 
% Las operaciones se ejecutan exactamente como en cualquier calculadora científica, declarando unos datos (en unos "contenedores" llamados variables) 
% y operando sobre ellos. Algunos ejemplos: 

% Scritps vs Command Window: FIGHT!
% Breakpoints y secciones. Use them!
% Variables: nuestros contenedores. Nombres y palabras reservadas.
% Errores. Aprovechemos el intérprete!
% Un espacio de trabajo limpio: clc y clear.

% Esto es un primer comentario. 

% Aritmética modular 
variable_1 = 1;                     %Primera variable de prueba ¿Qué pasa si no ponemos ;? %¿Por qué comentar el código y como hacerlo?
variable_2 = 2;                     %Segunda variable de prueba 
%Declara aquí una tercera variable
variable_3 = 10e3;

sum = variable_1+variable_2;       %Suma y guarda el resultado en sum
dif = variable_1-variable_2;       %Resta y guarda el resultado en diff
pro = variable_1*variable_2;       %Multiplica y guarda el resultado en pro 
div = variable_1/variable_2;       %Divide y guarda en div

%Con jerarquía de operaciones (utilizando paréntesis), calcula primero la suma de sum y dif con la tercera variable, 
%eleva al cuadrado el resultado y después hazle la raíz quinta.
result = ((sum+dif+variable_3)^2)^(1/5);
disp('Vamos a calcular el resultado!');
fprintf("El resultado de las anteriores operaciones es %.f \n", result);
fprintf("El resultado ha sido correcto");

%%
%Operaciones más avanzadas
pot = variable_1^variable_2;       %Potencias 
sq = sqrt(variable_1);             %Raíz cuadrada (solo la cuadrada!) 
variable_1 = -1;                   %Reasignación de valor a la variable variable_1
sqC = sqrt(variable_1);            %Raíz cuadrada (en el plano complejo!)

%%
%¿Cómo se calculan las funciones seno o coseno? ----> https://es.mathworks.com/

% Operaciones con matrices (the best!)
% Para MATLAB, todos los objetos que se manejan (prácticamente) son  tratados como arrays multidimensionales 
% (un escalar es una matriz de 1x1, un vector un array de nx1 o 1xn...). Esto hace que las operaciones matriciales 
% sean mucho más sencillas y rápidas.

%Definición de arrays
scalar = 2;
vector = [1; 1; 2];
matrix = [1 2 3 4; 8 0 0 0];

matrix_1 = rand(scalar, scalar);
matrix_2 = eye(scalar);
matrix_3 = magic(scalar);

%Indexación de arrays 
prueba = vector(1); 
prueba_2 = vector(1:2);
prueba_3 = matrix_1(:,1); 
prueba_4 = matrix_2(1,1);

concatenacion = [matrix_1 matrix_2; matrix_3 zeros(scalar, scalar)];

%Opera ahora
%ExProd(1:3,1) = scalar*vector;                 %Calcula el producto de un vector y un escalar
%ExProd(:,2) = ExProd(1,:)-vector;              %Algunas operaciones

%InProd(1) = dot(ExProd(1:2,2), vector(1:2,1));      %Producto escalar entre dos vectores
%InProd(2) = dot(ExProd(2,:), ExProd(2,:));          %Producto escalar del mismo vector consigo mismo (distancia euclídea!)
%InProd(3) = norm(ExProd(2,:))^2;                    %Producto escalar del mismo vector consigo mismo (distancia euclídea!)

MatrixProd = matrix_1*matrix_2*matrix_3;            %Producto de matrices
Error = matrix_1*vector(1:2,1);                     %Comprueba las dimensiones!

%%
%Cosas chulas de álgebra 
vector = vector.';                                  %Traspuesta de un vector/matriz
vector = vector.^scalar;                            %Operación vectorizada (componente a componente)
Det = det(matrix_1);                                %Determinante de una matriz
InvMatrix1 = inv(matrix_1);                         %Inversa de una matriz 
InvMatrix2 = matrix_1^(-1);                         %Inversa de una matriz. Y si ponemos .^? Observa la diferencia en el cómputo numérico!
[W, eigs] = eig(InvMatrix2);                        %Análisis del espectro de una matrix

%%
%Definiciones importantes 
xmax = 1;                                           %Longitud máxima
dx = 1e-3;                                          %Paso 
vector_1 = linspace(-5, 5, 100);                      %División regular en N pasos del intervalo -5 a 5
vector_2 = 0:dx:xmax;                               %División regular
vector_3 = [1 2 3];                                 %Vector

[X1, X2, X3] = meshgrid(vector_1, vector_2, vector_3);  %Producto cartesiano o mallado del espacio

% Algo de polinomios
pol = [1 2 3]; 
pol2 = poly(pol);
solutions = roots(pol2);

%%
% Ploteo y gráficas: https://www.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html
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
legend('Recta', 'Parábola'); 
xlabel('Eje x');
ylabel('Eje y');
title('Gráfica de prueba'); 
grid on;

%%
% Ejercicio 0
% Sergio tiene 3 gatos y dos perros. En total, pesan 30 kg. Los dos perros menos un gato
% pesan 10 kg. Todos los perros pesan lo mismo, igual para los gatos. Resuleve para ambos pesos en forma matricial.
A = [3 2; -1 2];
b = [30; 10];

% Si ahora dos gatos y cuatro perros pesan 7.5 kilos, calcula la solución de mínimos cuadrados.
A = [A; 2 4]; 
b = [b; 7.5]; 
sol = pinv(A)*b;

%%
% Ejercicio 1
% Representar la trayectoria de un avión que describe una hélice de radio 10, ascendiendo por el eje z.
% PISTA: para pasar 3 coordenadas, utiliza la función plot3 en vez de plot.
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
% potencial (aerodinámico, eléctrico, gravitatorio (?)). El espacio de configuración es de 2 x 2 metros (R2). 
% Conocemos la posición de tres fuentes (una carga positiva), de dos sumideros (carga negativa), en (cos(45º), -sin(45º))
% y un vórtice. Las cargas positivas son de 1,2 y 3 u; las negativas, -1, -3, -5; el vórtice posee vorticidad 10.

% Posición de las cargas: (0,1), (1,0), (0,-1)
% Posición de los sumideros: (cos(45º), -sin(45º)), (-1.5, -1)
% Posición de un vórtice: (0,0)

% Campo generado por una carga positiva: (qi/4*pi) * (1/||r-r_0||^2)
% Campo generado por una carga negativa: (qi/4*pi) * (1/||r-r_0||^2)
% Campo generado por un vórtice:    (gamma*i/2*pi) * (log z-z_0)

% Obtén la representación gráfica de las líneas de campo. 
L = 2;                      %Límite del espacio de configuracion
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

phi = real(F_t);                           %Definición de función potencial
psi = imag(F_t);                           %Definición de función de corriente

%Representación
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