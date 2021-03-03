%% Curso Crash Course de Matlab %% 
% Asociaci�n Aeroespacial Cosmos, URJC % 
% Madrid, febrero de 2021 % 

%% Introducci�n a MATrix LABoratory %%
% MATLAB 
% MATLAB es una abreviatura de MATrix LABoratory, un sistema de c�mputo num�rico que ofrece un entorno de desarrollo integrado 
% (IDE) con un lenguaje de programaci�n propio (lenguaje M).

% Si no tienes acceso a �l, y deseas utilizarlo sin conexi�n a internet, puedes descargar el programa con una licencia de prueba
% como estudiante: https://www.mathworks.com/campaigns/products/trials.html. Disponible para  Windows, macOS, Linux y Unix.

% EL ENTORNO
% Una vez abierto el programa, existen dos formas principales de comunicarse con el entorno: a trav�s de la ventana de comandos 
% (command window) o a trav�s de de la ventana del script. Cualquier informaci�n que generemos quedar� registrada en el workspace. La ventana
% Current Folder muestra la ra�z y los directorios en los que estamos trabajando.

% EL LENGUAGE
% Como decimos, el lenguaje propio de MATLAB se denomina lenguaje M. Este lenguaje es interpretado, din�mico, de alto nivel, y puede ejecutarse tanto 
% en el entorno interactivo, como a trav�s de un archivo de script (archivos *.m). Orientado al tratamiento de matrices y c�lculo num�rico,
% aunque se puede hacer pr�cticamente de todo. 

% Existen algunos entornos que trabajan con lenguaje M, con especial  menci�n a OCTAVE UPM y SCYLA, de licencia GNU. 
% JULIA o Python no son compatibles con M, pero son muy similares.

% TOOLBOXES Y DLC's
% Matlab est� construido por paquetes o toolboxes; cada uno de ellos contiene una serie de herramientas ya implementadas 
% especializadas en un campo: optimizaci�n matem�tica, tratamiento de se�ales, ingenier�a aeroespacial... La mayor�a de ellas 
% no est�n incluidas en los paquetes b�sicos, como el de estudiante. Especial relevancia tiene SIMULINK, una plataforma de an�lisis
% multidominio/mult�fisica de forma muy intuitiva

%% Parte 1: Matlab como calculadora 
% Algo de teor�a 
% Las operaciones se ejecutan exactamente como en cualquier calculadora cient�fica, declarando unos datos (en unos "contenedores" llamados variables) 
% y operando sobre ellos. Algunos ejemplos: 

% Scritps vs Command Window: FIGHT!
% Breakpoints y secciones. Use them!

% Aritm�tica modular 
variable_1 = ;                    %Primera variable de prueba �Qu� pasa si no ponemos ;? %�Por qu� comentar el c�digo y como hacerlo?
variable_2 = ;                     %Segunda variable de prueba 
%Declara aqu� una tercera variable

sum = variable_1+variable_2;       %Suma y guarda el resultado en sum
dif = variable_1-variable_2;       %Resta y guarda el resultado en diff
pro = variable_1*variable_2;       %Multiplica y guarda el resultado en pro 
div = variable_1/variable_2;       %Divide y guarda en div

%Con jerarqu�a de operaciones (utilizando par�ntesis), calcula primero la suma de sum y dif con la tercera variable, 
%eleva al cuadrado el resultado y despu�s hazle la ra�z quinta.
result = 
%fprintf("El resultado de las anteriores operaciones es %4.f \n", result);

%Operaciones m�s avanzadas
pot = variable_1^variable_2;       %Potencias 
sq = sqrt(variable_1);             %Ra�z cuadrada (solo la cuadrada!) 
variable_1 = -1;                   %Reasignaci�n de valor a la variable variable_1
sqC = sqrt(variable_1);            %Ra�z cuadrada (en el plano complejo!)

%�C�mo se calculan las funciones seno o coseno? ----> https://es.mathworks.com/

% Operaciones con matrices (the best!)
% Para MATLAB, todos los objetos que se manejan (pr�cticamente) son  tratados como arrays multidimensionales 
% (un escalar es una matriz de 1x1, un vector un array de nx1 o 1xn...). Esto hace que las operaciones matriciales 
% sean mucho m�s sencillas y r�pidas.

scalar = ;
vector = []; 
matrix_1 = rand(scalar, scalar);
matrix_2 = eye(scalar);
matrix_3 = magic(scalar);

%Preallocation
ExProd = zeros(size(vector,1), size(vector,2));     %Espera, espera. �Qu� prealocamos? �Qu� hace size()?
InProd = zeros(1,3);

%Opera ahora
ExProd(1,1:scalar) = scalar*vector;                 %Calcula el producto de un vector y un escalar
ExProd(2,:) = ExProd(1,:)-vector;                   %Algunas operaciones

InProd(1) = dot(ExProd(2,:), vector);               %Producto escalar entre dos vectores
InProd(2) = dot(ExProd(2,:), ExProd(2,:));          %Producto escalar del mismo vector consigo mismo (distancia eucl�dea!)
InProd(3) = norm(ExProd(2,:))^2;                    %Producto escalar del mismo vector consigo mismo (distancia eucl�dea!)

MatrixProd = matrix_1*matrix_2*matrix_3;            %Producto de matrices
Error = vector*matrix_1;                            %Comprueba las dimensiones!

%Cosas chulas de �lgebra 
vector = vector.';                                  %Traspuesta de un vector/matriz
vector = vector.^scalar;                            %Operaci�n vectorizada (componente a componente)
Det = det(matrix_1);                                %Determinante de una matriz
InvMatrix1 = inv(matrix_2);                         %Inversa de una matriz 
InvMatrix2 = matrix_2^(-1);                         %Inversa de una matriz. Y si ponemos .^?

%Definiciones importantes 
xmax = 1;                                           %Longitud m�xima
dx = 1e-3;                                          %Paso 
vector_1 = linspace();                              %Divisi�n regular
vector_2 = 0:dx:xmax;                               %Divisi�n regular
vector_3 = [];                                      %Un vector aleatorio

[X1, X2, X3] = meshgrid(vector_1, vector_2, vector_3);  %Producto cartesiano o mallado del espacio

% Ejercicio 1
% Sergio tiene 3 gatos y dos perros. En total, pesan 30 kg. Los dos perros menos un gato
% pesan 10 kg. Todos los perros pesan lo mismo, igual para los gatos. Resuleve para ambos pesos en forma matricial.


% Si ahora dos gatos pesan 7.5 kilos, calcula la soluci�n de m�nimos cuadrados.

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

%% Parte 2: Matlab como entorno de programaci�n 
% Algo de teor�a 

% Work out 

% Ejercicio 2

%% Parte 3: Matlab como calculadora avanzada 
% Algo de teor�a 

% Work out 

% Ejercicio 3

