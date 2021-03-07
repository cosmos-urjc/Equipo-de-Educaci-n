%% Parte 3: Matlab como calculadora avanzada 
% Algo de teor�a 
% Un campo vectorial es la distribuci�n de un vector en un cierto dominio,
% como la aceleraci�n gravitatoria en cada punto del espacio. Para
% nosotros, en general ser� un array de NxM dimensiones. 

% Matlab ofrece dos ventajas en el tratamiento de estos objetos, derivado de
% que puede trabajar f�cilmente con arrays: representaci�n y tratamiento.

% En muchas ocasiones, los campos vectoriales de estudio en ingenier�a y
% ciencia en general son funciones temporales, cuya evoluci�n se rige por
% ecuaciones diferenciales: las �rbitas terrestres, Navier-Stokes, etc... Es
% importante conocer algunas t�cnicas de integraci�n num�rica: Euler, RK4,
% odeXX (familia de solvers de Matlab)

% En ocasiones, habr� que integrar simult�neamente un conjunto de
% condiciones iniciales muy grande. Si utilizamos las herramientas
% proporcionadas por Matlab para el tratamiento de campos vectoriales,
% podemos reducir c�bicamente el tiempo que tarda el programa en integrar
% las trayectorias (en el espacio de estados). En general, ode45 no podr�
% utilizarse: solo admite campos vectoriales unidimensionales que propaga en
% el tiempo. Otros m�todos como RK4 ser�n necesarios

% Work out 

% Ejercicio: Ecuaci�n de Burguers 
L = 10;             %Espacio de configuraci�n (dominio, malla)
dx = 0.01;          %Paso espacial
x = -L:dx:L;        %Dominio de integraci�n

dt = 0.025;         %Paso temporal
tf = 100;           %Tiempo de integraci�n

nu = 0.01;          %Viscosidad del fluido, par�metro de la funci�n a integrar

%Condici�n inicial
f = sech(x);        %Condiciones iniciales

%Integraci�n y animaci�n
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
tf = 4000;              %Tiempo final de integraci�n
tspan = 0:dt:tf-dt;     %Span de integraci�n

%Par�metro gravitacional de la Tierra 
mu = 3.986e14;

%Condiciones iniciales
y0 = [7*10^6; 7*10^6; 0; 0; 2*10^3; 0];

%Acotamos el error de la integraci�n 
RelTol = 1e-14; 
AbsTol = 1e-14;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%Integramos con Ode45
[t, X] = ode45(@(t,y)kepler(t,y,mu), tspan, y0, options);

%Representaci�n
plot3(X(:,1),X(:,2),X(:,3));
grid on;

%Ejemplo de integraci�n vectorizada: Sistema de Lorenz
%Par�metros del sistema
beta = [10;8/3;28];

%Conjunto de condiciones iniciales (por ejemplo, part�culas fluidas)
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

%Integraci�n 
tf = 10; 
dt = 0.01; 
tspan = 0:dt:tf-dt; 

% Ejemplo 1: el sistema de Lorenz para un conjunto de condiciones iniciales
% sin c�digo vectorizado
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
    %Resoluci�n de la ecuaci�n de Burguers
    N = length(u); 
    fu = fft(u);                       %Transformada de Fourier del campo de velocidades
    k = (2*pi/L)*[-N/2:N/2-1];         %C�lculo de la longitud de onda espacial
    k = fftshift(k.');                 %Reordenamos la longitud de onda (arm�nicos)
    dfu = 1i*k.*fu;                    %Derivada espacial del campo de velocidades en el dominio de Fourier
    d2fu = -(k.^2).*fu;                %Segunda derivada espacial del campo de velocidades en el dominio de Fourier

    du = ifft(dfu);                    %Transformada inversa de la derivada del campo de velocidades
    d2u = ifft(d2fu);                  %Transformada inversa de la segunda derivada del campo de velocidades

    duout = -u.*du + nu*d2u;           %Ecuaci�n de Burguers con t�rmino advectivo y convectivo
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
    %Par�metros
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