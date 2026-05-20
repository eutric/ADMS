clear
close all
clc

f_max = 6;
eta = 1.5;
OM_max = f_max*2*pi;
lmax = @(EJ, m, eta, OM_max) pi*sqrt(sqrt(EJ/m)/eta/OM_max);
% Red
m_red = 122;
EA_red = 3.21e9;
EJ_red = 2e8;
lmax_red = lmax(EJ_red, m_red, eta, OM_max); % 14.9488
% Green
m_green = 80;
EA_green = 1.74e9;
EJ_green = 4.8e7;
lmax_green = lmax(EJ_green, m_green, eta, OM_max); % 11.6272
% Blue
m_blue = 400;
EA_blue = 2.5e9;
EJ_blue = 4.0e9;
lmax_blue = lmax(EJ_blue, m_blue, eta, OM_max); % 23.4930

% Geometrical data, A is the peak of first red beam, B of second one
% C is base of second red beam
% b is end of first green beam
coeff_bigblue = -(60-17.5)/67.6; % -.6287

% Known
R = 112;

% Red Beams points (upper case)
O = [0;  0];
A = [0;  60];
B = [17; A(2) + coeff_bigblue*17];        % [17; 49.3121]
C = [17; 0];
D = [50.60; 0];
E = [50.6; A(2) + coeff_bigblue*50.6];   % [50.6; 28.1879]
F = [67.6; 17.5];
G = [67.6; 0];
H = [128;0];
I = [67.6+60.4; R];

% Green beam points (lower case)
b = [C(1); B(2)-45];           % [17; 4.3121]
c = [O(1); 15];
d = [C(1); b(2) + c(2)];
e = [O(1); c(2) + c(2)];
f = [C(1); d(2) + c(2)];
g = [O(1); e(2) + c(2)];
h = [D(1); E(2)-15-8.75];
j = [G(1); F(2)/2];
k = [D(1); E(2)-15];

arco = @(theta) [
    I(1) + R*cos(theta);
    I(2) + R*sin(theta)
];
thetaF = acos((F(2)-I(2))/R) + pi/2;
thetaH = 3/2*pi;
arco_coord = arco(linspace(thetaF, thetaH, 4));

figure
% red 1
plot_beam(O, A, 'r')
hold on
plot_beam(C, B, 'r')
% red 2
plot_beam(D, E, 'r')
plot_beam(G, F, 'r')
% blue
plot_beam(A, F, 'b')
plot(arco_coord(1,:), arco_coord(2,:), 'b', LineWidth=1.5)
% green 1
plot_beam(O, b, 'g')
plot_beam(b, c, 'g')
plot_beam(c, d, 'g')
plot_beam(d, e, 'g')
plot_beam(e, f, 'g')
plot_beam(f, g, 'g')
plot_beam(g, B, 'g')
% green 2
plot_beam(G, h, 'g')
plot_beam(h, j, 'g')
plot_beam(j, k, 'g')
plot_beam(k, F, 'g')
grid on
n19 = mid_node(b,c);
scatter(n19(1), n19(2))
%% nodi della fem ('mazzati giulia)
n = [
    O(1), O(2);
    0, 7.5;
    0, 15;
    0, 22.5;
    0, 30;
    0, 37.5;
    0, 45;
    0, 52.5;
    0, 60;
    17, 0;
    17, 4.3121;
    17, 11.8121;
    17, 19.3121;
    17, 26.8121;
    17, 34.3121;
    17, 41.8121;
    17, 49.3121;
    mid_node(O,b)';
    mid_node(b,c)';
    mid_node(c,d)';
    mid_node(d,e)';
    mid_node(e,f)';
    mid_node(f,g)';
    mid_node(g,B)';
    mid_node(B,E,1/3)';
    mid_node(B,E,2/3)';
    D';
    h';
    k';
    mid_node(k,E)';
    E';
    G';
    j';
    mid_node(j,F)';
    F';
    mid_node(G,h)';
    mid_node(h,j)';
    mid_node(j,k)';
    mid_node(k,F)';
    arco_coord(:,2)';
    arco_coord(:,3)';
    arco_coord(:,4)';
    
];

norm(c-n19)

%% Struttura
clear
close all
clc

[xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure('Input.inp');
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
%% functions
function plot_beam(A, B, color)
    plot([A(1), B(1)], [A(2), B(2)], color, LineWidth=1.5);
end


