clearvars
close all

%% CLACULATION AND REPRESENTATION OF BEZIER CURVES
% Using two diferent methods, visualizing first the drawing of a curve
% using De Casteljau's algorithm, and then the simple calculation of the
% curve.

% Define number of steps desired
steps = 30;

% For animations, set to 1 (will save gif in folder
animation = 1;
imsize = 400;

% Define margins for graph
margin = 0.1;


%% POINT SELECTION POLYGON

Fdecas = figure(1); x0=200; y0=200; width=imsize; height=0.93*imsize;
set(Fdecas,'position',[x0,y0,width,height])
axis([0 1 0 1]);
title('Click to place points, press "Enter" when finished')
[x, y] = getpts();
clf(1,'reset')
P = [x'; y']';

%% Point checker
if size(P,2)<2
    error('At least two points must be defined')
elseif size(P,1)>size(P,2)
    P = P';
end

%% Preliminary definitions

% Order of curve
n = length(P)-1;

% Time vector
step = 1/steps;
tvec = linspace(0,1,steps);

%% DECASTELJAU ALGORITM ANIMATION

Fdecas = figure(1); x0=200; y0=200; width=imsize; height=0.93*imsize;
set(Fdecas,'position',[x0,y0,width,height])
hold on
title(sprintf('Bezier curve of order %i with De Casteljau algorithm', n));
init(P,margin);

Bdecas = B_casteljau(P,tvec);
im = cell(1,steps);

for i=1:steps
    t = tvec(i);
    
    % Delete the plotted temporary lines and dots
    if i>1
        delete(linemarker)
        delete(dotmarker)
    end
    
    % Plot the Bezier curve from 0 to the calculated instant
    plot(Bdecas(1,(1:i)),Bdecas(2,(1:i)),'k','LineWidth',3);
    
    % Plot the intermediate points
    [linemarker,dotmarker] = intermediate(P,t);
    
    if animation==1
        frame = getframe(Fdecas);
        im{i} = frame2im(frame);
    end
end

if animation==1
    writegif('decas.gif',im,steps,2/steps);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function initializes graph for Bezier curve.
% Plots given points and lines connectiong them
% Option to increase graph margins

function init(points,margin)

% Find edge points for axis scaling
M = max(points');
m = min(points');

% Define axes and margins
axis([(m(1)-margin) (M(1)+margin) (m(2)-margin) (M(2)+margin)]);

% Colour for polygon lines
c = gray;
c = c(30,:)/255;
c = [210, 211, 214]./255;

% Plot the initial points and their connecting lines
plot(points(1,:),points(2,:),'Color',c,'LineWidth',3);
scatter(points(1,:),points(2,:),'filled','k');
hold on

end

%% Function for calculation of Bezier curve 
% Using De Casteljau's algorithm

function curve = B_casteljau(points,tvec)

% Number of points
L = length(points);

% Find number of steps (time increments)
steps = length(tvec);

% Initialize variable for n-order curve
B = zeros(2,steps);

% Use De Casteljau's algorithm to calculate the curve at each interval
for i = 1:steps
    t = tvec(i);
    P = casteljau(points,L,t);
    B(:,i) = P(:,1,L);
end

curve = B;

end


%% This function uses De Casteljau's algorithm to 
% evaluate the set of points in each instant "t", and calculate the
% auxilairy points and lines
% The outputs to the function are the plotted data, which can be deleted
% to clear graph.

function [linemarker, dotmarker] = intermediate(points,t)

% Number of points
L = size(points,2);

% Extract points to the first iteration (defined points)
P(:,:,1) = points;

% Use time variable as marker
xlabel(sprintf('t = %.2f', t));

% Generate colour gradient for plots
c = parula(L);

% Calulate Bezier curve (recursive)
P = casteljau(points,L,t);

% Plotting in a double iteration
for i=1:L
    for j=1:L-i
        
        % Plot auxiliary points
        if i==(L-1)
            dotmarker(i,j) = scatter(P(1,j,i+1),P(2,j,i+1),40,c(i,:),'filled');
        else
            dotmarker(i,j) = scatter(P(1,j,i+1),P(2,j,i+1),10,c(i,:));
        end

        % Plot the lines between points
        if (j > 1) & (L > 2)
            linemarker(i,j) = plot(P(1,j-1:j,i+1),P(2,j-1:j,i+1),'Color',c(i,:));
        elseif (L < 3)
            linemarker = [];
        end
        
    end
end

end

%% Function to write the animation frames to a gif file

function writegif(filename, im, nImages, delay)

for i = 1:nImages
    [A,map] = rgb2ind(im{i},256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end

end


%% Iterative portion of De Casteljau's algorithm

function P = casteljau(points,L,t)

% Extract points to the first iteration (defined points)
P(:,:,1) = points;

for i=1:1:L
    for j=1:1:L-i
        % Calulate Bezier curve (recursive)
        P(:,j,i+1) = (1-t)*P(:,j,i) + t*P(:,j+1,i); 
    end
end

end