clear
close all
clc

% Cantilever Beam
syms A B C D x omega gamma E J L

w_x =  A*cos(gamma*x) + B*sin(gamma*x) + C*cosh(gamma*x) + D*sinh(gamma*x);
dw_x = diff(w_x, x)
ddw_x = diff(dw_x, x)
dddw_x = diff(ddw_x, x)

% BCs
% w_x(0) = 0, w'_x(0) = 0, EJ*w''_x(L) = 0, EJ*w'''_x(L) = 0

H = [
    1, 0, 1, 0;
    0, gamma, 0, gamma;
    -E*J*gamma^2*cos(gamma*L), -E*J*gamma^2*sin(gamma*L), E*J*gamma^2*cosh(gamma*L), E*J*gamma^2*sinh(gamma*L);
    E*J*gamma^3*sin(gamma*L), -E*J*gamma^3*cos(gamma*L), E*J*gamma^3*sinh(gamma*L), E*J*gamma^3*cosh(gamma*L);
];

detH = matlabFunction(det(H));

% Dati
L = 1200e-3; % m
h = 8e-3; %    m
b = 40e-3; %   m
rho = 2700; %  kg/m^3
E = 68e9; %    Pa

A = b*h;
m = rho*A;

J = b*h^3/12;

gamma = @(omega) (m*omega.^2./(E*J)).^(1/4);

fmax = 500; % Hz

detH = @(omega) detH(E,J,L,gamma(omega));

figure
fs = linspace(0, fmax, 1000*fmax);
omegas = 2*pi*fs;
gammas = gamma(omegas);
abs_detH = abs(detH(omegas));
semilogy(fs, abs_detH)
grid on

neg_detH = -abs_detH;
[nat_freq, locs] = findpeaks(neg_detH);
fs(locs)

hold on

scatter(fs(locs),abs_detH(locs))