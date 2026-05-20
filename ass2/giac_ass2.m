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
lmax_blue =lmax(EJ_blue, m_blue, eta, OM_max); % 23.4930

% Geometrical data, A is the peak of first red beam, B of second one
% C is base of second red beam
% b is end of first green beam (potrà esse il modo più comodo???)
coeff_bigblue = -(60-17.5)/67.6; % -.6287

% Known
R = 112;
% Red Beams point (upper case)
O.x = 0;
O.y = 0;
A.x = 0;
A.y = 60;
B.x = 17;
B.y = A.y + coeff_bigblue*B.x; % 49.3121
C.x = 17;
C.y = 0;
D.x = 50.60;
D.y = 0;
E.x = 50.6;
E.y = A.y + coeff_bigblue*E.x; % 28.1879
F.x = 67.6;
F.y = 17.5;
G.x = 67.6;
G.y = 0;
I.x = 67.6+60.4;
I.y = R;

% green beam points, lower case
b.x = C.x;
b.y = B.y-45; % 4.3121
c.x = O.x;
c.y = 15;
d.x = C.x;
d.y = b.y + c.y;
e.x = O.x;
e.y = c.y+c.y;
f.x = C.x;
f.y = d.y + c.y;
g.x = O.x;
g.y = e.y+c.y;
h.x = D.x;
h.y = E.y-15-8.75;
j.x = G.x;
j.y = F.y/2;
k.x = D.x;
k.y = E.y-15;



arco = @(theta) [
    I.x + R*cos(theta);
    I.y + R*sin(theta)
];

thetaF = acos((F.y-I.y)/R)+pi/2;
thetaH = 3/2*pi;

arco_coord = arco(linspace(thetaF,thetaH, 100));

figure
% red 1
plot_beam(O,A,'r')
hold on
plot_beam(C,B,'r')
% red 2
plot_beam(D,E,'r')
plot_beam(G,F,'r')
% blue
plot_beam(A,F,'b')
plot(arco_coord(1,:),arco_coord(2,:),'b', LineWidth=1.5)
% green 1
plot_beam(O,b,'g')
plot_beam(b,c,'g')
plot_beam(c,d,'g')
plot_beam(d,e,'g')
plot_beam(e,f,'g')
plot_beam(f,g,'g')
plot_beam(g,B,'g')
% green 2
plot_beam(G,h,'g')
plot_beam(h,j,'g')
plot_beam(j,k,'g')
plot_beam(k,F,'g')

grid on

%% nodi della fem ('mazzati giulia)
b2 = sqrt(b.x^2+b.y^2)
btheta = atan(b.y/b.x);
n18.x = b2/2*cos(btheta);
n18.y = b2/2*sin(btheta); 
bc = sqrt((c.x-b.x)^2+(c.y-b.y)^2)
cd = sqrt((d.x-c.x)^2+(d.y-c.y)^2)
%% functions
function plot_beam(A,B, color)
    plot([A.x, B.x], [A.y, B.y], color, lineWidth = 1.5);
end