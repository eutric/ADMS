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
clear 
close all
clc

% Beam Data

beam.L=1.200;  % m
beam.h=0.008;  % m
beam.b=0.040;  % m
beam.rho=2700; % kg/m^3
beam.E=68e9;   % Pa - Alluminio
beam.m=beam.rho*beam.h*beam.b; % kg/m

beam.J=beam.h^3*beam.b/12;     % m^4;

% Freq data
f_max = 200;  % Hz
f_res = 10000; % Hz
f_vect = linspace(0, f_max, f_res);
om_vect = 2*pi.*f_vect; % rad/s
g = @(om) (beam.m*om.^2/beam.E/beam.J).^(1/4);

%build the H matrix
H=@(g) [1 0 1 0;
    0 g 0 g;
    -beam.E*beam.J.*g.^2.*cos(g*beam.L) -beam.E*beam.J.*g.^2.*sin(g*beam.L) beam.E*beam.J.*g.^2.*cosh(g*beam.L) beam.E*beam.J.*g.^2.*sinh(g*beam.L);
    beam.E*beam.J.*g.^3.*sin(g*beam.L) -beam.E*beam.J*g.^3.*cos(g*beam.L) beam.E*beam.J*g.^3.*sinh(g*beam.L) beam.E*beam.J*g.^3.*cosh(g*beam.L) ];
dets=zeros(length(om_vect),1);
for i = 1:length(om_vect)
    dets(i)=det(H(g(om_vect(i))));
end
abs_det = abs(dets);
[~, i_nat] = findpeaks(-abs(dets));

figure
semilogy(f_vect, abs_det)
hold on
scatter (f_vect(i_nat), abs_det(i_nat))
grid on
% computing mode shapes
figure
xx=linspace(0,beam.L,1000);
labels = cell(1, length(i_nat));
for i =1:length(i_nat)
    mode(i).i_nat=i_nat(i);
    mode(i).OM=om_vect(i_nat(i));
    mode(i).G=g(mode(i).OM);
    mode(i).HH=H(mode(i).G);
    mode(i).f=mode(i).OM/2/pi;
    mode(i).X=[1;-mode(i).HH(2:end,2:end)\mode(i).HH(2:end,1)];
    mode(i).phi=@(x)mode(i).X(1)*cos(mode(i).G.*x) + mode(i).X(2)*sin(mode(i).G.*x)+ mode(i).X(3)*cosh(mode(i).G.*x) + mode(i).X(4)*sinh(mode(i).G.*x);

    plot(xx, mode(i).phi(xx))
    tag = sprintf('Mode %d', i);
    labels{i} = {tag};
    hold on
end
grid on
legend([labels{:}])

% compute FRF
%assume 
x_in=.2; % m
x_out=.6;

[G,mode]=FRF_num(mode,x_in,x_out,beam,-1); % L'ultimo parametro è la fase/contro fase di spostamento e forza, la prof li ha di segno opposto
figure
subplot(2,1,1)
semilogy(f_vect,abs(G(om_vect)),LineWidth=2)
grid on
subplot(2,1,2)
plot(f_vect,angle(G(om_vect)),LineWidth=2)
grid on

% test FRF in other point, compute FRF_num
x_out_vect=[0.2,0.7,1.2];
[FRF1,mode1]=FRF_num(mode,x_in,x_out_vect(1),beam);
[FRF2,mode2]=FRF_num(mode,x_in,x_out_vect(2),beam);
[FRF3,mode3]=FRF_num(mode,x_in,x_out_vect(3),beam);
FRF_list={FRF1,FRF2,FRF3};


% work backwords compute parameters starting from sperimental FRFs
%assume aproximation around mode 1

% x = [om_i,csi_1,A1,A2,A3,RL1,RL2,RL3,RH1,RH2,RH3];

range = 1; %Hz

[x_sols, G_num_list, x_0s, n] = mode_i_num(FRF_list, range, f_vect, mode);

G_in_out = cell(n.ch,1);
for jj=1:n.ch
    G_in_out{jj} = @(OM) 0;
    for ii =1:n.modes
        G_in_out{jj} = @(OM) G_in_out{jj}(OM) + G_num_list{ii,jj}(x_sols(:,ii), OM);
        f_vect_test(ii,:) = linspace(mode(ii).f-range, mode(ii).f+range, 50);
    end
end

om_vect_test = 2*pi*f_vect_test;

% A figure for each channel, near each mode
for jj=1:n.ch
    fig_name = sprintf('Canale %d', jj);
    figure('Name', fig_name, 'NumberTitle', 'off');
    for ii=1:n.modes
        subplot(2,n.modes,ii)
        semilogy(f_vect_test(ii,:), abs(FRF_list{jj}(om_vect_test(ii,:))))
        hold on
        semilogy(f_vect_test(ii,:), abs(G_in_out{jj}(om_vect_test(ii,:))),'o')
        grid on
        tit = sprintf('Modo %d - f = %.2f Hz', ii, x_0s(1,ii)/2/pi);
        title(tit)
        subplot(2,n.modes,ii+n.modes)
        plot(f_vect_test(ii,:), angle(FRF_list{jj}(om_vect_test(ii,:))))
        hold on
        plot(f_vect_test(ii,:), angle(G_in_out{jj}(om_vect_test(ii,:))), 'o')
        grid on
    end

end

% comparison between frf (exp - num)

figure
subplot(2,1,1)
semilogy(f_vect, abs(G_in_out{1}(om_vect)), 'o')
grid on
hold on
semilogy(f_vect, abs(FRF1(om_vect)));

subplot(2,1,2)
plot(f_vect, angle(G_in_out{1}(om_vect)), 'o')
grid on
hold on
plot(f_vect, angle(FRF1(om_vect)));

% comparison between mode shapes
% i need the first 4 numeric eigen vector

% point 5

% From fitting, i can, now, estimate the values of the mode shape and the 
% points where i evaluated the response

% 1. From mode_i_num(), i can now estimate A = phi_ij*phi_ik

% 2. Fixing the input (k), i can at first derive phi_ik for each mode and
%    then, for any output, derive phi_ij

% 3. The modeshape is then known for each output point: phi_ij

% 1 - x_j = x_k 
x_j = .65; % Location of the input

[FRF_jj{1}] = FRF_num(mode, x_j, x_j, beam, 1);
[x_sols, ~, ~, ~] = mode_i_num(FRF_jj, 1, f_vect, mode);

% phi_ij = zeros(n.modes,1);
Ajj = x_sols(3,:);
phi_ij = sqrt(real(Ajj));
% for ii = 1:n.modes
%     phi_ij(ii) = sqrt( 1i*2*x_sols(2,ii)*x_sols(1,ii)^2*G_num_list{ii,1}(x_sols(:,ii),x_sols(1,ii)) );
% end

n_outs = 20;
x_ks = linspace(0,beam.L, n_outs); % 6 output points, equally spaced

FRFs = cell(n_outs, 1);
for ii = 1:n_outs
    FRFs{ii} = FRF_num(mode, x_j, x_ks(ii), beam, -1);
end

[x_sols, G_num_list, x_0s, n] = mode_i_num(FRFs, 1, f_vect, mode);

% For each solution, i have the estimated A, from wich i can deriva 
% phi_ij
Ajk = x_sols(3:2+n.ch,:);
phi_ik = zeros(n.ch, n.modes);

for ii=1:n.modes
    phi_ik(:,ii) = real(Ajk(:,ii))./phi_ij(ii);
end

% for ii=1:n.modes
%     for jj=1:n.ch
%         phi_ik(ii,jj) = 1i*G_num_list{ii,jj}(x_sols(:,ii),x_sols(1,ii))*2*x_sols(2,ii)*x_sols(1,ii)^2/phi_ij(ii);
%     end
% end

% => A is matrix n.ch x n.modes, the element jj,ii is the i-th mode value
% on point j of the i-th modeshape
% => to reconstruct the shape of each mode, it's enough to plot these
% points (?)

figure
for ii=1:n.modes
    subplot(2,2,ii)
    plot(xx, mode(ii).phi(xx));
    hold on
    if (mode(ii).phi(x_ks(end))>0 && phi_ik(end,ii) > 0) || (mode(ii).phi(x_ks(end))<0 && phi_ik(end,ii) < 0)
        scatter(x_ks, phi_ik(:,ii))
    else
        scatter(x_ks, -phi_ik(:,ii))
    end
    grid on
end

%% what happens around .6?
clear
close all
clc

% Beam Data

beam.L=1.200;  % m
beam.h=0.008;  % m
beam.b=0.040;  % m
beam.rho=2700; % kg/m^3
beam.E=68e9;   % Pa - Alluminio
beam.m=beam.rho*beam.h*beam.b; % kg/m

beam.J=beam.h^3*beam.b/12;     % m^4;

% Freq data
f_max = 200;  % Hz
f_res = 10000; % Hz
f_vect = linspace(0, f_max, f_res);
om_vect = 2*pi.*f_vect; % rad/s
g = @(om) (beam.m*om.^2/beam.E/beam.J).^(1/4);

%build the H matrix
H=@(g) [1 0 1 0;
    0 g 0 g;
    -beam.E*beam.J.*g.^2.*cos(g*beam.L) -beam.E*beam.J.*g.^2.*sin(g*beam.L) beam.E*beam.J.*g.^2.*cosh(g*beam.L) beam.E*beam.J.*g.^2.*sinh(g*beam.L);
    beam.E*beam.J.*g.^3.*sin(g*beam.L) -beam.E*beam.J*g.^3.*cos(g*beam.L) beam.E*beam.J*g.^3.*sinh(g*beam.L) beam.E*beam.J*g.^3.*cosh(g*beam.L) ];
dets=zeros(length(om_vect),1);
for i = 1:length(om_vect)
    dets(i)=det(H(g(om_vect(i))));
end
abs_det = abs(dets);
[amp, i_nat] = findpeaks(-abs(dets));

labels = cell(1, length(i_nat));
for i =1:length(i_nat)
    mode(i).i_nat=i_nat(i);
    mode(i).OM=om_vect(i_nat(i));
    mode(i).G=g(mode(i).OM);
    mode(i).HH=H(mode(i).G);
    mode(i).f=mode(i).OM/2/pi;
    mode(i).X=[1;-inv(mode(i).HH(2:end,2:end))*mode(i).HH(2:end,1)];
    mode(i).phi=@(x)mode(i).X(1)*cos(mode(i).G.*x) + mode(i).X(2)*sin(mode(i).G.*x)+ mode(i).X(3)*cosh(mode(i).G.*x) + mode(i).X(4)*sinh(mode(i).G.*x);
end


x_in = .8;
x_out = .3429;

[FRF, mode] = FRF_num(mode, x_in, x_out, beam, -1);

figure
subplot(2,1,1)
semilogy(f_vect, abs(FRF(om_vect)))
grid
subplot(2,1,2)
plot(f_vect, angle(FRF(om_vect)))
grid