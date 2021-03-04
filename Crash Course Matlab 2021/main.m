%% Curso Crash Course de Matlab %% 
% Asociación Aeroespacial Cosmos, URJC % 
% Madrid, febrero de 2021 % 

%% Introducción a MATrix LABoratory %%
% MATLAB 
% MATLAB es una abreviatura de MATrix LABoratory, un sistema de cómputo numérico que ofrece un entorno de desarrollo integrado 
% (IDE) con un lenguaje de programación propio (lenguaje M).

% Si no tienes acceso a él, y deseas utilizarlo sin conexión a internet, puedes descargar el programa con una licencia de prueba
% como estudiante: https://www.mathworks.com/campaigns/products/trials.html. Disponible para  Windows, macOS, Linux y Unix.

% EL ENTORNO
% Una vez abierto el programa, existen dos formas principales de comunicarse con el entorno: a través de la ventana de comandos 
% (command window) o a través de de la ventana del script. Cualquier información que generemos quedará registrada en el workspace. La ventana
% Current Folder muestra la raíz y los directorios en los que estamos trabajando.

% EL LENGUAGE
% Como decimos, el lenguaje propio de MATLAB se denomina lenguaje M. Este lenguaje es interpretado, dinámico, de alto nivel, y puede ejecutarse tanto 
% en el entorno interactivo, como a través de un archivo de script (archivos *.m). Orientado al tratamiento de matrices y cálculo numérico,
% aunque se puede hacer prácticamente de todo. 

% Existen algunos entornos que trabajan con lenguaje M, con especial  mención a OCTAVE UPM y SCYLA, de licencia GNU. 
% JULIA o Python no son compatibles con M, pero son muy similares.

% TOOLBOXES Y DLC's
% Matlab está construido por paquetes o toolboxes; cada uno de ellos contiene una serie de herramientas ya implementadas 
% especializadas en un campo: optimización matemática, tratamiento de señales, ingeniería aeroespacial... La mayoría de ellas 
% no están incluidas en los paquetes básicos, como el de estudiante. Especial relevancia tiene SIMULINK, una plataforma de análisis
% multidominio/multífisica de forma muy intuitiva

%% Parte 1: Matlab como calculadora 
% Algo de teoría 
% Las operaciones se ejecutan exactamente como en cualquier calculadora científica, declarando unos datos (en unos "contenedores" llamados variables) 
% y operando sobre ellos. Algunos ejemplos: 

% Scritps vs Command Window: FIGHT!
% Breakpoints y secciones. Use them!
% Variables: nuestros contenedores. Nombres y palabras reservadas.
% Errores. Aprovechemos el intérprete!

% Aritmética modular 
variable_1 = ;                     %Primera variable de prueba ¿Qué pasa si no ponemos ;? %¿Por qué comentar el código y como hacerlo?
variable_2 = ;                     %Segunda variable de prueba 
%Declara aquí una tercera variable

sum = variable_1+variable_2;       %Suma y guarda el resultado en sum
dif = variable_1-variable_2;       %Resta y guarda el resultado en diff
pro = variable_1*variable_2;       %Multiplica y guarda el resultado en pro 
div = variable_1/variable_2;       %Divide y guarda en div

%Con jerarquía de operaciones (utilizando paréntesis), calcula primero la suma de sum y dif con la tercera variable, 
%eleva al cuadrado el resultado y después hazle la raíz quinta.
result = 
disp('Vamos a calcular el resultado!');
%fprintf("El resultado de las anteriores operaciones es %4.f \n", result);

%Operaciones más avanzadas
pot = variable_1^variable_2;       %Potencias 
sq = sqrt(variable_1);             %Raíz cuadrada (solo la cuadrada!) 
variable_1 = -1;                   %Reasignación de valor a la variable variable_1
sqC = sqrt(variable_1);            %Raíz cuadrada (en el plano complejo!)

%¿Cómo se calculan las funciones seno o coseno? ----> https://es.mathworks.com/

% Operaciones con matrices (the best!)
% Para MATLAB, todos los objetos que se manejan (prácticamente) son  tratados como arrays multidimensionales 
% (un escalar es una matriz de 1x1, un vector un array de nx1 o 1xn...). Esto hace que las operaciones matriciales 
% sean mucho más sencillas y rápidas.

%Definición de arrays
scalar = ;
vector = []; 
matrix_1 = rand(scalar, scalar);
matrix_2 = eye(scalar);
matrix_3 = magic(scalar);

%Indexación de arrays 
prueba = vector(1); 
prueba_2 = vector(1:2);
prueba_3 = matrix_1(2,:); 
prueba_4 = matrix_2(1,1);

concatenacion = [matrix_1 matrix_2; matrix_3 eye(5)];

%Preallocation
ExProd = zeros(size(vector,1), size(vector,2));     %Espera, espera. ¿Qué prealocamos? ¿Qué hace size()?
InProd = zeros(1,3);

%Opera ahora
ExProd(1,1:scalar) = scalar*vector;                 %Calcula el producto de un vector y un escalar
ExProd(2,:) = ExProd(1,:)-vector;                   %Algunas operaciones

InProd(1) = dot(ExProd(2,:), vector);               %Producto escalar entre dos vectores
InProd(2) = dot(ExProd(2,:), ExProd(2,:));          %Producto escalar del mismo vector consigo mismo (distancia euclídea!)
InProd(3) = norm(ExProd(2,:))^2;                    %Producto escalar del mismo vector consigo mismo (distancia euclídea!)

MatrixProd = matrix_1*matrix_2*matrix_3;            %Producto de matrices
Error = vector*matrix_1;                            %Comprueba las dimensiones!

%Cosas chulas de álgebra 
vector = vector.';                                  %Traspuesta de un vector/matriz
vector = vector.^scalar;                            %Operación vectorizada (componente a componente)
Det = det(matrix_1);                                %Determinante de una matriz
InvMatrix1 = inv(matrix_2);                         %Inversa de una matriz 
InvMatrix2 = matrix_2^(-1);                         %Inversa de una matriz. Y si ponemos .^? Observa la diferencia en el cómputo numérico!
[W, eigs] = eig(InvMatrix2);                        %Análisis del espectro de una matrix

%Definiciones importantes 
xmax = 1;                                           %Longitud máxima
dx = 1e-3;                                          %Paso 
vector_1 = linspace(-5, 5, N);                      %División regular en N pasos del intervalo -5 a 5
vector_2 = 0:dx:xmax;                               %División regular
vector_3 = [];                                      %Un vector aleatorio

[X1, X2, X3] = meshgrid(vector_1, vector_2, vector_3);  %Producto cartesiano o mallado del espacio

% Algo de polinomios
pol = []; 
pol2 = poly(MatrixProd);
solutions = roots(pol2);

% Ploteo y gráficas: https://www.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html
figure(i) 
hold on 
hold off 
legend(); 
xlabel();
title(); 
view(3);

% Ejercicio 0
% Sergio tiene 3 gatos y dos perros. En total, pesan 30 kg. Los dos perros menos un gato
% pesan 10 kg. Todos los perros pesan lo mismo, igual para los gatos. Resuleve para ambos pesos en forma matricial.


% Si ahora dos gatos pesan 7.5 kilos, calcula la solución de mínimos cuadrados.

% Ejercicio 1
% Representar la trayectoria de un avión que describe una hélice de radio 10, ascendiendo por el eje z.
% PISTA: para pasar 3 coordenadas, utiliza la función plot3 en vez de plot.
% z = 0:pi/50:10*pi;
% x = sin(t);               
% y = cos(t);
% plot3(x,y,z)

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
x = -L:h:L-h;               %Dominio horizontal
y = -L:h:L-h;               %Dominio vertical
[X, Y] = meshgrid(x,y);     %Mallado
Z = X-1i*Y;                 %Isomorfismo entre R2 y C

F_v = ((1i*Gamma)/(2*pi))*(log(z-b*1i));   %Potencial torbellino
F_f = ((Q)/(2*pi))*(log(z+a-b*1i));        %Potencial fuente
F_t = F_v+F_f;

phi = real(F_t);                           %Definición de función potencial
psi = imag(F_t);                           %Definición de función de corriente

%Representación
figure(1)
plot3(X,Y,E_horizontal);  

figure(2) 
dens = 20;
streamslice(X,Y,E_horizontal,E_vertical,dens);

figure(3) 
hold on
contour(X,Y,E_horizontal,100);
contour(X,Y,E_vertical,100);
hold off

figure(4) 
Z = sin(X);
surf(X,Y,Z);

%% Parte 2: Matlab como entorno de programación 
% Algo de teoría: ¿Qué ofrece Matlab como lenguaje de programación? 

% Work out 
% Más palabras reservadas!
% Tipo boolean. Comparaciones. Funciones any, all, indexación booleana
% Bucles: for, while. 
% Sentencias condicionadas: 
% Estructura switch! 
% Estructuras y listas de cells
% Funciones. ¿Qué son y para qué sirven? ¿Cómo se implementan?


% Ejercicio 2
% Se desea implementar la primera versión software de un controlador de
% ángulo de ataque para un perfil sustentador NACA experimental (Cl_alpha = 2*pi/sqrt(1-M^2), modelo compresible)  
% que irá implementado en una aeronave también experimental de carga y 
% transporte de ayuda humanitaria. Se desea en todo momento, durante el
% vuelo de crucero y a medida que la carga se va desplegando, mantener una
% configuración de vuelo horizontal rectílineo uniforme.

% Durante esta fase del vuelo, el controlador de empuje corrige siempre las
% variaciones de resistencia, de manera que el problema longitudinal no es
% relevante. La condición, por tanto, a asegurar será 
% eqn(1) --> L = W

% Debido a que se trata de una aeronave experimental, el controlador solo
% podrá generar incrementos de ángulo de ataque aleatorios uniformemente distribuidos 
% en el intervalo de 0 a delta\alpha, donde delta\alpha 
% es el cambio de ángulo necesario para compensar el peso.  

% Asumáse una polar sustentadora solo dependiente del ángulo de ataque (CL_0 = 0). 
% Él ángulo de ataque inicial es de 2 deg (asumáse que el perfil no entra en pérdida nunca).
% La densidad y aceleración gravitatoria son constantes, a altitud de 5000m. La velocidad de crucero
% de diseño es Mach 0.65 y es constante. La superfice alar son 25 m^2. 
% La atmósfera se asume regida por modelo ISA. Finalmente, cada 5
% segundos, la aeronave pierde 50 kg. La duración total de la fase de
% crucero es de 300s y el peso en despegue es de 5 toneladas. 

% Para corroborar el diseño aerodinámico, es necesario que como output del
% programa se muestre el ángulo de ataque en función del tiempo.

% Datos 
%Variables independientes
h = 5000; 
g = 9.81;
tf = 300; 
dt = 5; 
t = 0:dt:tf-dt; %Vector de tiempo

%Características de la aeronave
S = 25; 
Mach = 0.65;
M0 = 5*10^3; 
deltaM = -50;

%Características del perfil 
Cl_0 = 0; 
Cl_alpha = 2*pi/sqrt(1-Mach^2); %Teoría potencial linealizada; analogía de Prandt-Glauert
alpha0 = 2*pi/180;

% Resolución 
%Step 1: cálculo la densidad y temperatura con el modelo ISA 
rho = 0.736116;
T = 255.650;
Rg = 287; 
gamma = 1.4;

%Step 2: cálcula la velocidad de vuelo y la presión dinámica en la fase de
%crucero 
V = Mach*sqrt(Rg*gamma*T);
p_d = (1/2)*rho*V^2;
alpha(1) = alpha0;

%Step 3: cálculo principal en estructura de bucle temporal
for i = 1:size(t,2)-1
    W(i) = g*(M0+deltaM*i);
    Cl = Cl_0+Cl_alpha*alpha(i);
    alpha_n(i) = (W(i))/(p_d*S*Cl);
    if (alpha_n(i) == 0)
        alpha(i+1) = alpha(i);
    elseif (alpha_n(i)-alpha(i) > 0)&&(alpha_n(i) ~= 0)
        alpha(i+1) = alpha(i) + abs(alpha_n(i)-alpha(i))*rand;
    elseif (alpha_n(i)-alpha(i) <= 0)&&(alpha_n(i) ~= 0)
        alpha(i+1) = alpha(i) - abs(alpha_n(i)-alpha(i))*rand;
    end
end

% Plotting 
%Representa el ángulo de ataque como función del tiempo 
figure(2)
plot(t,alpha); 
xlabel('Tiempo (s)')
ylabel('Ángulo de ataque (rad)')

% Animación del perfil rotando
%Se desea realizar una animación del giro del perfil con el tiempo (como
%cabecea o encabrita). Para ello, para cada instante de tiempo, se debe
%definir una matriz de rotación Q que gire todo el contorno del perfil un
%ángulo alpha 

%Contorno del perfil 
c = 1;       %Cuerda del perfil adimensionalizada
e  = 0.13;   %Espesor en tanto por uno
dx = 0.0001; %Paso geométrico de la malla
x =0:dx:c; 	 %Dominio de representación del perfil 
    
ye = 5*e*(0.2969.*sqrt(x)-0.126.*(x)-0.3515.*(x.^2)+0.2843*(x.^3)-0.1015*(x.^4));    %Ley de espesor
yc = (2429241163689825*((4.*x)/3 - 1).*(((4.*x)/3 - 1).^2 - 1))/72057594037927936;   %Ley de curvatura
yti = yc + ye; %Forma final de intradós
yte = yc - ye; %Forma final de extradós

%Bucle principal 
figure(1)
for i = 1:size(t,2)
    clf
    Q(1:3,1:3) = [cos(alpha(i)) sin(alpha(i)) 0; -sin(alpha(i)) cos(alpha(i)) 0; 0 0 1];
    for j = 1:size(x,2)
        V_aux1 = Q*[x(j);yte(j);0];
        V_aux2 = Q*[x(j);yti(j);0];
        Extrados_girado(j) = V_aux1(2);
        Intrados_girado(j) = V_aux2(2);
    end
    hold on 
    plot(x, Extrados_girado,'b');
    plot(x, Intrados_girado,'b');
    hold off
    axis([0 1 -1 1]);
    xlabel ('Componente x del perfil')
    ylabel ('Componente y del perfil')
    pause(0.7);
end

%% Parte 3: Matlab como calculadora avanzada 
% Algo de teoría 

% Work out 

% Ejercicio 3

