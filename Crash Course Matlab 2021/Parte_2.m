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