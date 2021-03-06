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
% Un espacio de trabajo limpio: clc y clear.

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
a < b 
b > a 
a => b 
b <= a
a == b 
a != b

true && false 
false && true 
true && true 
true || false

vector1 = randi([0 100], 1, 20);
vector2 = rand(1,20);
comparacion = (vector1 == vector2);
comparacion2 = (any(vector1) == 1); 
comparacion3 = (all(vector2) <= 0);
comparacion4 = (max(vector1) >= max(vector2)); 

a = true; 
b = false; 

if (a)
elseif (b)
else
end

if (a > b)
else
end

a = 'Hola'; 
b = 'Adios'; 

switch (a)
    case ('Hola')
    case ('Adios')
    otherwise 
end

% Bucles: for, while. 
num = 100; 
n = 0;
for i = 1:num
    n = n+1;
end

for i = 1:num
    A(i,:) = 0;
    for j = i:num 
        A(i,j) = 1;
    end
end

t = 0:1e-1:10;
for i = 1:num 
    y(i,:) = i*sin(t);
end

while (a < 1)
end

GoOn = true;
i = 0;
while (true) && (GoOn)
    i = i+1;
    if (i > 3)
        GoOn = false;
    end
end

% Estructuras y listas de cells
Estructura.Edad = [10 20 30 40]; 
Estructura.Nombres = ['Carlos' 'Juan' 'Alicia']; 
Estructura.Media = sum(Estructura.Edad)/size(Estructura.Edad); 
fprintf('La media de edad es de %.4f años \n', Estructura.Media); 

lista = {}; 
lista{1} = 'Carlos'; 
lista{2} = 10; 
lista{3} = true; 

% Funciones. ¿Qué son y para qué sirven? ¿Cómo se implementan?
Eqn1 = @(x)(x^2+2*x+1);
Matrix = @(x)([1 2 3 x; 4 5 6 7]);
procedimiento();
[Media] = calcula_media(Estructura.Edad);

% Ejercicio 3
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
% Un campo vectorial es la distribución de un vector en un cierto dominio,
% como la aceleración gravitatoria en cada punto del espacio. Para
% nosotros, en general será un array de NxM dimensiones. 

% Matlab ofrece dos ventajas en el tratamiento de estos objetos, derivado de
% que puede trabajar fácilmente con arrays: representación y tratamiento.

% En muchas ocasiones, los campos vectoriales de estudio en ingeniería y
% ciencia en general son funciones temporales, cuya evolución se rige por
% ecuaciones diferenciales: las órbitas terrestres, Navier-Stokes, etc... Es
% importante conocer algunas técnicas de integración numérica: Euler, RK4,
% odeXX (familia de solvers de Matlab)

% En ocasiones, habrá que integrar simultáneamente un conjunto de
% condiciones iniciales muy grande. Si utilizamos las herramientas
% proporcionadas por Matlab para el tratamiento de campos vectoriales,
% podemos reducir cúbicamente el tiempo que tarda el programa en integrar
% las trayectorias (en el espacio de estados). En general, ode45 no podrá
% utilizarse: solo admite campos vectoriales unidimensionales que propaga en
% el tiempo. Otros métodos como RK4 serán necesarios

% Work out 

% Ejercicio: Ecuación de Burguers 
L = 10;             %Espacio de configuración (dominio, malla)
dx = 0.01;          %Paso espacial
x = -L:dx:L;        %Dominio de integración

dt = 0.025;         %Paso temporal
tf = 100;           %Tiempo de integración

nu = 0.01;          %Viscosidad del fluido, parámetro de la función a integrar

%Condición inicial
f = sech(x);        %Condiciones iniciales

%Integración y animación
for i = 1:tf
    elaps = i*dt; 
    [~,u] = ode45(@(t,u)burguers(t,u,L,nu), [0 dt], f);
    f = u(end,:);
    plot(x, real(f))
    axis([-L L -1.5 1.5])
    pause(0.01)
end

% Ejercicio 3: Ecuaciones de Kepler  
dt = 0.1;               %Paso temporal
tf = 4000;              %Tiempo final de integración
tspan = 0:dt:tf-dt;     %Span de integración

%Parámetro gravitacional de la Tierra 
mu = 3.986e14;

%Condiciones iniciales
y0 = [7*10^6; 7*10^6; 0; 0; 2*10^3; 0];

%Acotamos el error de la integración 
RelTol = 1e-14; 
AbsTol = 1e-14;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%Integramos con Ode45
[t, X] = ode45(@(t,y)kepler(t,y,mu), tspan, y0, options);

%Representación
plot3(X(:,1),X(:,2),X(:,3));
grid on;

%Ejemplo de integración vectorizada: Sistema de Lorenz
%Parámetros del sistema
beta = [10;8/3;28];

%Conjunto de condiciones iniciales (por ejemplo, partículas fluidas)
L = 10;
dx = 1; 
dy = 1; 
dz = 1;

xvec = -L:dx:L; 
yvec = -L:dy:L; 
zvec = -L:dz:L; 
[x0,y0,z0] = meshgrid(xvec,yvec,zvec);

X(1,:,:,:) = x0; 
X(2,:,:,:) = y0; 
X(3,:,:,:) = z0; 

plot3(X(1,:),X(2,:),X(3,:), 'r.', 'LineWidth', 4, 'MarkerSize',4);
axis([-40 40 -40 40 -40 40]);

%Integración 
tf = 10; 
dt = 0.01; 
tspan = 0:dt:tf-dt; 

% Ejemplo 1: el sistema de Lorenz para un conjunto de condiciones iniciales
% sin código vectorizado
% for i = 1:size(tspan,2)
%     t = dt*i; 
%     for j = 1:size(xvec,2) 
%         for k = 1:size(yvec,2)
%             for z = 1:size(zvec,2)
%                 y_in = X(:,j,k,z); 
%                 y_out = rk4singlestep(@(t,y)loren_dynamics(t,y,beta), dt, t, y_in);
%                 y_trajec(:,j,k,z) = y_out;
%             end
%         end
%     end
%     plot3(y_trajec(1,:), y_trajec(2,:), y_trajec(3,:), 'r.', 'LineWidth', 2, 'MarkerSize',4);
%     axis([-40 40 -40 40 -40 40]);
% end

y_in = X;
for i = 1:size(tspan,2)
    t = dt*i; 
    y_trajec = rk4singlestep(@(t,y)lorenz_data(t,y,beta), dt, t, y_in); 
    y_in = y_trajec; 
    plot3(y_trajec(1,:), y_trajec(2,:), y_trajec(3,:), 'r.', 'LineWidth', 1, 'MarkerSize',1);
    axis([-40 40 -40 40 -40 40]);
    view(20,40)
    drawnow
end

%% Funciones auxiliares 
function [duout] = burguers(~, u, L, nu)
    %Resolución de la ecuación de Burguers
    N = length(u); 
    fu = fft(u);                       %Transformada de Fourier del campo de velocidades
    k = (2*pi/L)*[-N/2:N/2-1];         %Cálculo de la longitud de onda espacial
    k = fftshift(k.');                 %Reordenamos la longitud de onda (armónicos)
    dfu = 1i*k.*fu;                    %Derivada espacial del campo de velocidades en el dominio de Fourier
    d2fu = -(k.^2).*fu;                %Segunda derivada espacial del campo de velocidades en el dominio de Fourier

    du = ifft(dfu);                    %Transformada inversa de la derivada del campo de velocidades
    d2u = ifft(d2fu);                  %Transformada inversa de la segunda derivada del campo de velocidades

    duout = -u.*du + nu*d2u;           %Ecuación de Burguers con término advectivo y convectivo
end

function [dy] = kepler(mu, t, X)
    %Variables a integrar como coordenadas de nuestro vector X   
    r = X(1:3).';
    v = X(4:6).';
    R = norm(r);
    
    %Ecuaciones del movimiento kepleriano
    gamma = -mu*(r/R^3);
    dy = [v gamma].'; 
end

function [dy] = lorenz_data(t,y,beta)
    %Parámetros
    sigma = beta(1);
    gamma = beta(2); 
    rho =  beta(3); 

    %Variables 
    X = y(1,:,:,:); 
    Y = y(2,:,:,:); 
    Z = y(3,:,:,:);

    dy = [ sigma*(Y-X); 
           X.*(rho-Z)-Y; 
           X.*Y-gamma*Z];
end

function [y_out] = rk4singlestep(fun, dt, initial_time, initial_position)
    f1 = fun(initial_time, initial_position); 
    f2 = fun(initial_time + dt/2, initial_position+(dt/2)*f1); 
    f3 = fun(initial_time + dt/2, initial_position+(dt/2)*f2); 
    f4 = fun(initial_time + dt,   initial_position+dt*f3);
    y_out = initial_position + (dt/6)*(f1+2*f2+2*f3+f4);
end