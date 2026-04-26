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


%%
%countilever beam parameters
clear all
close all
clc
beam.L=1.200;  % mm
beam.h=0.008;
beam.b=0.040;
beam.rho=2700; %kg/mm^3
beam.E=68e9; %MPa
beam.m=beam.rho*beam.h*beam.b; % kg/mm
beam.J=beam.h^3*beam.b/12; % mm^4;
f_max = 200;
f_res = 10000;
f_vect=linspace(0,f_max,f_res);
om_vect=2*pi.*f_vect;
g=@(om)(beam.m*om.^2/beam.E/beam.J).^(1/4);

%build the H matrix
H=@(g)[1 0 1 0;
    0 g 0 g;
    -beam.E*beam.J.*g.^2.*cos(g*beam.L) -beam.E*beam.J.*g.^2.*sin(g*beam.L) beam.E*beam.J.*g.^2.*cosh(g*beam.L) beam.E*beam.J.*g.^2.*sinh(g*beam.L);
    beam.E*beam.J.*g.^3.*sin(g*beam.L) -beam.E*beam.J*g.^3.*cos(g*beam.L) beam.E*beam.J*g.^3.*sinh(g*beam.L) beam.E*beam.J*g.^3.*cosh(g*beam.L) ];
dets=[];
for i = 1:length(om_vect)
    dets(i)=det(H(g(om_vect(i))));
end
semilogy(om_vect,abs(dets))
hold on
i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
        
    end
end
semilogy(om_vect(i_nat),abs(dets(i_nat)),'or')
grid on
%natural frequencies
%computing mode shapes
figure
xx=linspace(0,beam.L,1000);
for i =1:length(i_nat)
    mode(i).i_nat=i_nat(i);
    mode(i).OM=om_vect(i_nat(i));
    mode(i).G=g(mode(i).OM);
    mode(i).HH=H(mode(i).G);
    mode(i).f=mode(i).OM/2/pi;
    mode(i).X=[1;-inv(mode(i).HH(2:end,2:end))*mode(i).HH(2:end,1)];
    mode(i).phi=@(x)mode(i).X(1)*cos(mode(i).G.*x) + mode(i).X(2)*sin(mode(i).G.*x)+ mode(i).X(3)*cosh(mode(i).G.*x) + mode(i).X(4)*sinh(mode(i).G.*x);
    plot(xx,mode(i).phi(xx))
    hold on

end
grid on
%% compute FRF
%assume 
x_in=.2; % m
x_out=1.2;

[G,mode]=FRF_num(mode,x_in,x_out,beam,-1); % L'ultimo parametro è la fase/contro fase di spostamento e forza, la prof li ha di segno opposto
figure
subplot(2,1,1)
semilogy(f_vect,abs(G(om_vect)),LineWidth=2)
grid on
subplot(2,1,2)
plot(f_vect,angle(G(om_vect)),LineWidth=2)
grid on
    
%% test FRF in other point, compute FRF_num
x_out_vect=[0.2,0.5,0.7];
[FRF1,mode1]=FRF_num(mode,x_in,x_out_vect(1),beam);
[FRF2,mode2]=FRF_num(mode,x_in,x_out_vect(2),beam);
[FRF3,mode3]=FRF_num(mode,x_in,x_out_vect(3),beam);
FRF_list={FRF1,FRF2,FRF3};
figure
subplot(2,3,1)
semilogy(f_vect,abs(FRF1(om_vect)),LineWidth=2);
grid
subplot(2,3,2)
semilogy(f_vect,abs(FRF2(om_vect)),LineWidth=2);
grid
subplot(2,3,3)
semilogy(f_vect,abs(FRF3(om_vect)),LineWidth=2);
grid
subplot(2,3,4)
plot(f_vect,angle(FRF1(om_vect)),LineWidth=2);
grid
subplot(2,3,5)
plot(f_vect,angle(FRF2(om_vect)),LineWidth=2);
grid

subplot(2,3,6)
plot(f_vect,angle(FRF3(om_vect)),LineWidth=2);
grid

%% work backwords compute parameters starting from sperimental FRFs
%assume aproximation around mode 1

G_exp = [FRF1(om_vect)',FRF2(om_vect)', FRF3(om_vect)'];
% x = [om_i,csi_1,A1,A2,A3,RL1,RL2,RL3,RH1,RH2,RH3];

G_num = @(x,OM) x(3)/(-OM.^2+1i*2*x(2)*x(1).*OM+x(1)^2)+x(6)/OM^2+x(9);
G_num1 =   @(x,OM) x(4)/(-OM.^2+1i*2*x(2)*x(1).*OM+x(1)^2)+x(7)/OM^2+x(10);
G_num2 =   @(x,OM) x(5)/(-OM.^2+1i*2*x(2)*x(1).*OM+x(1)^2)+x(8)/OM^2+x(11);


range = 1; %Hz
for ii = 1:3
[mode_ampl(:,ii),mode_loc_exp(:,ii)] = findpeaks(abs(G_exp(:,ii)));
end


err =@(x) [G_exp(mode_loc_exp(1)-f_res/f_max * range:mode_loc_exp(1)+f_res/f_max * range,1) - G_num1(x,f_vect(mode_loc_exp(1)-f_res/f_max * range:mode_loc_exp(1)+f_res/f_max * range));
    G_exp(mode_loc_exp(2)-f_res/f_max * range:mode_loc_exp(2)+f_res/f_max * range,2) - G_num1(x,f_vect(mode_loc_exp(2)-f_res/f_max * range:mode_loc_exp(2)+f_res/f_max * range));
    G_exp(mode_loc_exp(3)-f_res/f_max * range:mode_loc_exp(3)+f_res/f_max * range,3) - G_num1(x,f_vect(mode_loc_exp(3)-f_res/f_max * range:mode_loc_exp(3)+f_res/f_max * range));   
];


% x


x01 = [mode_loc_exp(1,1),];
tol = 1e-4;
k = find(abs(abs(G_exp(:,1)) - mode_ampl(1,1)/sqrt(2)) < tol);





