%% Parte 2: Matlab como entorno de programaci�n 
% Algo de teor�a: �Qu� ofrece Matlab como lenguaje de programaci�n? 

% Work out 
% M�s palabras reservadas!
% Tipo boolean. Comparaciones. Funciones any, all, indexaci�n booleana
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
fprintf('La media de edad es de %.4f a�os \n', Estructura.Media); 

lista = {}; 
lista{1} = 'Carlos'; 
lista{2} = 10; 
lista{3} = true; 

% Funciones. �Qu� son y para qu� sirven? �C�mo se implementan?
Eqn1 = @(x)(x^2+2*x+1);
Matrix = @(x)([1 2 3 x; 4 5 6 7]);
procedimiento();
[Media] = calcula_media(Estructura.Edad);

% Ejercicio 3
% Se desea implementar la primera versi�n software de un controlador de
% �ngulo de ataque para un perfil sustentador NACA experimental (Cl_alpha = 2*pi/sqrt(1-M^2), modelo compresible)  
% que ir� implementado en una aeronave tambi�n experimental de carga y 
% transporte de ayuda humanitaria. Se desea en todo momento, durante el
% vuelo de crucero y a medida que la carga se va desplegando, mantener una
% configuraci�n de vuelo horizontal rect�lineo uniforme.

% Durante esta fase del vuelo, el controlador de empuje corrige siempre las
% variaciones de resistencia, de manera que el problema longitudinal no es
% relevante. La condici�n, por tanto, a asegurar ser� 
% eqn(1) --> L = W

% Debido a que se trata de una aeronave experimental, el controlador solo
% podr� generar incrementos de �ngulo de ataque aleatorios uniformemente distribuidos 
% en el intervalo de 0 a delta\alpha, donde delta\alpha 
% es el cambio de �ngulo necesario para compensar el peso.  

% Asum�se una polar sustentadora solo dependiente del �ngulo de ataque (CL_0 = 0). 
% �l �ngulo de ataque inicial es de 2 deg (asum�se que el perfil no entra en p�rdida nunca).
% La densidad y aceleraci�n gravitatoria son constantes, a altitud de 5000m. La velocidad de crucero
% de dise�o es Mach 0.65 y es constante. La superfice alar son 25 m^2. 
% La atm�sfera se asume regida por modelo ISA. Finalmente, cada 5
% segundos, la aeronave pierde 50 kg. La duraci�n total de la fase de
% crucero es de 300s y el peso en despegue es de 5 toneladas. 

% Para corroborar el dise�o aerodin�mico, es necesario que como output del
% programa se muestre el �ngulo de ataque en funci�n del tiempo.

% Datos 
%Variables independientes
h = 5000; 
g = 9.81;
tf = 300; 
dt = 5; 
t = 0:dt:tf-dt; %Vector de tiempo

%Caracter�sticas de la aeronave
S = 25; 
Mach = 0.65;
M0 = 5*10^3; 
deltaM = -50;

%Caracter�sticas del perfil 
Cl_0 = 0; 
Cl_alpha = 2*pi/sqrt(1-Mach^2); %Teor�a potencial linealizada; analog�a de Prandt-Glauert
alpha0 = 2*pi/180;

% Resoluci�n 
%Step 1: c�lculo la densidad y temperatura con el modelo ISA 
rho = 0.736116;
T = 255.650;
Rg = 287; 
gamma = 1.4;

%Step 2: c�lcula la velocidad de vuelo y la presi�n din�mica en la fase de
%crucero 
V = Mach*sqrt(Rg*gamma*T);
p_d = (1/2)*rho*V^2;
alpha(1) = alpha0;

%Step 3: c�lculo principal en estructura de bucle temporal
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
%Representa el �ngulo de ataque como funci�n del tiempo 
figure(2)
plot(t,alpha); 
xlabel('Tiempo (s)')
ylabel('�ngulo de ataque (rad)')

% Animaci�n del perfil rotando
%Se desea realizar una animaci�n del giro del perfil con el tiempo (como
%cabecea o encabrita). Para ello, para cada instante de tiempo, se debe
%definir una matriz de rotaci�n Q que gire todo el contorno del perfil un
%�ngulo alpha 

%Contorno del perfil 
c = 1;       %Cuerda del perfil adimensionalizada
e  = 0.13;   %Espesor en tanto por uno
dx = 0.0001; %Paso geom�trico de la malla
x =0:dx:c; 	 %Dominio de representaci�n del perfil 
    
ye = 5*e*(0.2969.*sqrt(x)-0.126.*(x)-0.3515.*(x.^2)+0.2843*(x.^3)-0.1015*(x.^4));    %Ley de espesor
yc = (2429241163689825*((4.*x)/3 - 1).*(((4.*x)/3 - 1).^2 - 1))/72057594037927936;   %Ley de curvatura
yti = yc + ye; %Forma final de intrad�s
yte = yc - ye; %Forma final de extrad�s

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