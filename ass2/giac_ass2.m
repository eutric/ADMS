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

%% FEM
clear
close all
clc

[xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure('Input.inp');
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);

% Adding spring
k_B = 5e7; 
n_B = 42;
% n_B = 18;
K(idb(n_B, 2), idb(n_B, 2)) = K(idb(n_B, 2), idb(n_B, 2)) + k_B;

% Considering structural damping
alpha = .1;   % 1/s
beta  = 2e-4; % s

C = alpha*M + beta*K;

% Modes
Mff = M(1:ndof, 1:ndof);
Mfc = M(1:ndof, ndof+1:end);
Mcf = Mfc';
Mcc = M(ndof+1:end, ndof+1:end);

Kff = K(1:ndof, 1:ndof);
Kfc = K(1:ndof, ndof+1:end);
Kcf = Kfc';
Kcc = K(ndof+1:end, ndof+1:end);

Cff = C(1:ndof, 1:ndof);
Cfc = C(1:ndof, ndof+1:end);
Ccf = Cfc';
Ccc = C(ndof+1:end, ndof+1:end);

[modes, om2] = eig(Mff\Kff); % Eigenvalues problem
om_i = sqrt(diag(om2)); % Natural freqs
[om_i,index] = sort(om_i);
f_i = om_i/2/pi;
modes = modes(:,index);

n_m = 4;
scale_factor = 20;

figure
for ii=1:n_m
    subplot(2,2,ii)
    diseg2(modes(:,ii),scale_factor,incidenze,l,gamma,posiz,idb,xy)
    tit = sprintf('Mode %d - f = %.4f Hz', ii, f_i(ii));
    title(tit)
end

% Forces
% n_A = 1; % free node
n_A = 9;
F0 = zeros(ndof,1);
F0(idb(n_A,1)) = 1; % Horizontal force on A

FRF = @(OM, F0) (-OM.^2*Mff + 1i*OM*Cff + Kff)\F0;
f_res = .01;
fs = 0:f_res:6;
OMs = 2*pi*fs;

% Computing FRFs for each freq
FRFs = zeros(ndof, length(OMs));
for ii = 1:length(OMs)
    FRFs(:,ii) = FRF(OMs(ii),F0);
end

% Requested OUTPUTs
FRF_to_xA = FRFs(idb(n_A,1),:);
FRF_to_yB = FRFs(idb(n_B,2),:);

figure
subplot(2,2,1)
semilogy(fs, abs(FRF_to_xA),'r',LineWidth=1.5)
hold on
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency Hz');
ylabel('| |');
title('Magnitude - Horizontal displacement on node A')
grid on
subplot(2,2,3)
plot(fs, angle(FRF_to_xA),'r',LineWidth=1.5)
hold on
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency Hz');
title('Phase - Horizontal displacement on node A')
grid on

subplot(2,2,2)
semilogy(fs, abs(FRF_to_yB),'r',LineWidth=1.5)
hold on
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency Hz');
ylabel('| |');
title('Magnitude - Vertical displacement on node B')
grid on
subplot(2,2,4)
plot(fs, angle(FRF_to_yB),'r',LineWidth=1.5)
hold on
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency Hz');
title('Phase - Vertical displacement on node B')
grid on


% 4 
% a - modal superposition
tic
n_m = 2; % I consider only 2 first modes
Phi = modes(:,1:n_m);

F0 = zeros(ndof, 1);
F0(idb(n_B,2)) = 1; % Vertical force on B

Mmod = Phi'*Mff*Phi;
Kmod = Phi'*Kff*Phi;
Cmod = Phi'*Cff*Phi;
Q0 = Phi'*F0; % Lagrangian component
q0 = @(OM,Q0) (-OM.^2*Mmod + 1i*OM*Cmod + Kmod)\Q0;

% Psi(:,1) = Phi(:,1)/sqrt(Mmod(1,1));
% Psi(:,2) = Phi(:,2)/sqrt(Mmod(2,2));
% Mmod_norm = Psi'*Mff*Psi; % experiment, not needed
% Kmod_norm = Psi'*Kff*Psi;
% Cmod_norm = Psi'*Cff*Psi;
% Q0_norm = Psi'*F0;
% q0 = @(OM,Q0) (-OM.^2*Mmod_norm + 1i*OM*Cmod_norm + Kmod_norm)\Q0_norm;


X0_modal = zeros(ndof, length(OMs));
q0s = zeros(n_m, length(OMs));
for ii = 1:length(OMs)
    q0s(:,ii) = q0(OMs(ii),Q0);
    X0_modal(:,ii) = Phi*q0s(:,ii);
end
toc
% b - FEM
tic
X0_4 = zeros(ndof, length(OMs));
for ii = 1:length(OMs)
    X0_4(:,ii) = FRF(OMs(ii), F0); 
end
toc
% Comparing
figure
subplot(1,2,1)
semilogy(fs, abs(X0_4(idb(n_A,1),:)),'r',LineWidth=1.5)
hold on
semilogy(fs, abs(X0_modal(idb(n_A,1),:)), 'm', LineWidth=1.5)
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency [Hz]');
ylabel('| |');
legend('FEM','I and II Mode superposition')
title('Magnitude - From Vert in B to Horiz in A')
grid on
subplot(1,2,2)
plot(fs, angle(X0_4(idb(n_A,1),:)),'r',LineWidth=1.5)
hold on
plot(fs, angle(X0_modal(idb(n_A,1),:)),'m',LineWidth=1.5)
xline([f_i(1), f_i(2), f_i(3), f_i(4)], 'g')
xlabel('Frequency [Hz]');
ylabel('Phase')
legend('FEM','I and II Mode superposition')
title('Phase - From Vert in B to Horiz in A')
grid on
% Va bene perché aumentando i modi, diventano uguali
% Da capire cosa succede a circa .91 Hz, aggiungendo i modi ad alta
% frequenza lì si crea un antirisonanza... non c'ha un vers, non so se è
% spiegabile però

% figure
% for ii=1:n_m
%     semilogy(fs, abs(q0s(ii,:)),LineWidth=1.5)
%     hold on
% end
% grid on

% 5 - Static deflection of a vertical distributed load
pG = [  % N/m - global reference
    0;
    -7000;
    0;
];    

el_blues = 23:30; % Number of elements ! NOT NODES

F = zeros(ndof,1);

Lambda = @(g) [ 
    cos(g), -sin(g), 0;
    sin(g),  cos(g), 0;
    0,       0,      1;
];

Lambda2 = @(g) [ 
    Lambda(g), zeros(3);
    zeros(3),  Lambda(g)
];

for k = el_blues
    pL = Lambda(gamma(k))'*pG;
    Fnodes = [
        l(k) * pL(1); 
        l(k)/2 * pL(2); 
        l(k)^2/12 * pL(2); 
        l(k) * pL(1)/2; 
        l(k)/2 * pL(2); 
        -l(k)^2/12 * pL(2)
    ];
    
    Fnodes_global = Lambda2(gamma(k))*Fnodes;
    
    E = zeros(length(Fnodes),ndof);
    E(:,incidenze(k,:)) = eye(length(Fnodes));
    
    F = F + E'*Fnodes_global;
end

x =  Kff\F;
figure()
diseg2(x,scale_factor,incidenze,l,gamma,posiz,idb,xy)


max_w = 0;
labels = {};
% figure
for k = el_blues 

    xglobal = x(incidenze(k,:));
    xlocal = Lambda2(gamma(k))'*xglobal;
    csi_vect= linspace(0,l(k),100000);
    % Cubic shape function for w(csi) - vertical displacement
    a = xlocal(2); 
    b = xlocal(3);
    c = -3/l(k)^2 * xlocal(2) + 3/l(k)^2 * xlocal(4) -2/l(k)*xlocal(3) -1/l(k)*xlocal(6); 
    d = 2/l(k)^3 * xlocal(2) -2/l(k)^3 *xlocal(4) +1/l(k)^2*xlocal(3)+1/l(k)^2*xlocal(6);
    
    coeff = [d,c,b,a];
    w = polyval(coeff,csi_vect);
    
    
    % plot(csi_vect,w, 'DisplayName', sprintf('Elemento %d', k));
    % labels{end+1} = sprintf('Elemento %d', k);
    % legend(labels); 
    % hold on
    
    [max_old,max_index] = max(abs(w));
    
    if max_w < max_old
        max_w = max_old;
        elemento = k;
        maxloc = max_index;
    end
end

v0 = 2; % m/s
a1 = 3.5;
a2 = 1.5;
cum = [0, cumsum(l(23:30))];

% Prima legge oraria, cerco t in cui: cum(27-21) = v01*t + 1/2*a1*t^2
v1 = @(t) v0 + a1*t;
a = 1/2*a1;
b = v0;
c = -cum(27-21);
t27 = (-b+sqrt(b^2-4*a*c))/(2*a);
v27 = v1(t27);

v2 = @(t) v27 + a2*t;
a = 1/2*a2;
b = v27;
c = -(cum(30-21)-cum(27-21));
t30 = t27 + (-b+sqrt(b^2-4*a*c))/(2*a);

x_t = @(t) (t<=t27)*(v0*t + 1/2*a1*t^2) + (t>t27 & t<=t30)*(cum(27-21)+v27*(t-t27)+1/2*a2*(t-t27)^2);

m = 90; % kg
F = [
    0;
    m*9.81; % Col segno più ha più senso - forza locale di là
    0;
];
n = length(M);
n_m = 10;
Phi = modes(:,1:n_m);

Mq = Phi'*Mff*Phi;
Kq = Phi'*Kff*Phi;
Cq = Phi'*Cff*Phi;

Mq_inv = Mq\eye(n_m);

A = [
    zeros(n_m), eye(n_m);
    -Mq_inv*Kq, -Mq_inv*Cq
];

Ft = @(t) [
    zeros(n_m,1);
    zeros(n_m,1)+(t<=t30)*Mq_inv*Phi'*F_t(F, x_t, t, n, ndof, l, gamma, incidenze)
];
odefun = @(t, z) A*z + Ft(t);
ts = linspace(0,50,1000);
[ts, zs] = ode45(odefun, ts, zeros(length(A),1));

figure
hold on
labels = {};
for ii=1:n_m
    plot(ts, zs(:,ii))
    labels{end+1} = sprintf('Modo %d', ii);
    legend(labels); 
end

xs = Phi*zs(:,1:n_m)';
figure
plot(ts, xs(idb(n_B,2),:))
grid on
for ii=1:length(ts)
    Fts(:,ii) = F_t(F, x_t, ts(ii), n, ndof, l, gamma, incidenze);
end
figure
plot(ts, Fts(idb(n_B,2),:))
%% functions

function plot_beam(A, B, color)
    plot([A(1), B(1)], [A(2), B(2)], color, LineWidth=1.5);
end
