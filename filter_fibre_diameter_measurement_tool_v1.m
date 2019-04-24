% Tool to help measure HEPA fibre diameters manually on-screen
%
% v1 24/04/19
% J B McClean
% Imperial College London
% j.mcclean15@imperial.ac.uk

clear
clc

%% Detect edges

% Read in the image (must be 1024 x 768)
I = imread('YOUR_IMAGE_HERE.tif');

% Generate a binary image of the edges
I_edges = edge(I,'Canny');

% Plot the results
I_both = imfuse(I,I_edges);

imtool(I);

%% Generate random points

imtool close all

% Define crosshair size
crosshair_length = 50;
crosshair_width = 2;
crosshair_offset = (crosshair_width - 1)/2;

% Abbreviate points to make it quicker to define vertices of crosshair
o = crosshair_offset;
l = crosshair_offset + crosshair_length;

% Define crosshair centre point
x_ch = randi(1024,1);
y_ch = randi(768,1);
crosshair_centre = [x_ch y_ch];

% Define extreme points of crosshair
crosshair_E = [x_ch + crosshair_offset + crosshair_length, y_ch];
crosshair_N = [x_ch, y_ch + crosshair_offset + crosshair_length];
crosshair_W = [x_ch - crosshair_offset - crosshair_length, y_ch];
crosshair_S = [x_ch, y_ch - crosshair_offset - crosshair_length];

% Calculate verticles of crosshair polygon
P_ch = [x_ch + o, y_ch + o;
        x_ch + o, y_ch + l;
        x_ch - o, y_ch + l;
        x_ch - o, y_ch + o;
        x_ch - l, y_ch + o;
        x_ch - l, y_ch - o;
        x_ch - o, y_ch - o;
        x_ch - o, y_ch - l;
        x_ch + o, y_ch - l;
        x_ch + o, y_ch - o;
        x_ch + l, y_ch - o;
        x_ch + l, y_ch + o];
    
P_x = P_ch(:,1);
P_y = P_ch(:,2);
I_crosshair = poly2mask(P_x, P_y, 768, 1024);

% Generate corners of viewport for zoomed-in view
view_x1 = x_ch - 2*crosshair_length;
view_x2 = x_ch + 2*crosshair_length;
view_y2 = y_ch + 2*crosshair_length;
view_y1 = y_ch - 2*crosshair_length;

% Check limits
if(view_x1 < 1)
    shift = -view_x1;
    view_x1 = 1;
    view_x2 = view_x2 + shift;
end

if(view_x2 > 1024)
    shift = view_x2 - 1024;
    view_x2 = 1024;
    view_x1 = view_x1 - shift;
end

if(view_y1 < 1)
    shift = -view_y1;
    view_y1 = 1;
    view_y2 = view_y2 + shift;
end

if(view_y2 > 768)
    shift = view_y2 - 768;
    view_y2 = 768;
    view_y1 = view_y1 - shift;
end

I_with_crosshair = imfuse(I, I_crosshair, 'method', 'blend');
hFig = imtool(I_with_crosshair(view_y1:view_y2, view_x1:view_x2), ...
    'InitialMagnification', 475);
set(hFig, 'units','normalized','outerposition',[0.05 0 0.95 1])

%% Read in diameters in pixels

d_f = [];

for i = 1:100
   disp(['d_f = [d_f distance', int2str(i), ']']);     
   eval(['d_f = [d_f distance', int2str(i), '];']); 
end

%% Covert diameters in pixels to diameters in microns

% You need to work this out

d_f_micron4 = d_f./(19 + 1/3);

%% Plot CDF of measured fibres for microglass

figure(4)
histogram(d_f_micron,100,'Normalization','cdf','DisplayStyle','stairs');
xlabel('Fibre diameter (micron)');
ylabel('Fraction of fibres');
ylim([0 1]);
xlim([0.1 2.3647]);
title('CDF of measured fibre diameter - central microglass layer');
hold on
% y_fit lognormal parameters values are:
% 3428 aa (image _02) = -0.3006, 0.8952
% 3248 aa (image _03) = -0.5189, 0.6730
% 3428 aa (image _04) = -0.7593, 0.7734
% images _02, _03 and _04 combined = -0.5262, 0.8053
y_fit = logncdf(linspace(0,8,1000),-0.5262,0.8053);
plot(linspace(0,8,1000),y_fit);
set(gca,'XScale','log');
grid on;

% Object to help with exporting data to .csv
[N, edges] = histcounts(d_f_micron,100,'Normalization','cdf');
figure(40)
stairs(edges(1:end-1),N);
set(gca, 'XScale', 'log');
xlim([0.1 2.3647]);
ylim([0 1]);

fileName = join(['moxie_3428aa_03_d_f.csv'])
EXPORTABLE = [edges(1:end-1)' N'];
csvwrite(fileName, EXPORTABLE);

%% Plot CDF of measured fibres for PET

figure(4)
histogram(d_f_micron,100,'Normalization','cdf','DisplayStyle','stairs');
xlabel('Fibre diameter (micron)');
ylabel('Fraction of fibres');
title('CDF of measured fibre diameter - exterior PET layer');
hold on
xlim([20 max(d_f_micron)]);
ylim([0 1]);
y_fit = logncdf(linspace(0,113,1000),3.5990,0.2456);
plot(linspace(0,113,1000),y_fit);
set(gca,'XScale','log');
grid on;

%% Fit a lognormal distribution

% Plot histogram of measured fibres with a fitted distribution
figure(2)
histfit(d_f_micron,20,'lognormal');
title('Histogram of measured fibre distribution');

% Get fitted distribution parameters
pd = fitdist(d_f_micron','lognormal');
mu = pd.mu;
sigma = pd.sigma;

% Plot a histogram of simulated fibres from the fitted distrib. fn.
fibre_diameters_simulated = [];
for i = 1:1000
   fibre_diameters_simulated(i) = lognrnd(mu,sigma);    
end
figure(3)
histfit(fibre_diameters_simulated,100,'lognormal');
title('Histogram of simulated fibre distribution');

mu_actual = exp(mu + 0.5*sigma^2);
sigma_actual = sqrt((exp(sigma^2) - 1)*(exp(2*mu + sigma^2)));





