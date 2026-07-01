clear
close all
clc
%load data
FRF_data = load("FRF_H1.mat");
f_vect = FRF_data.f;
FRFs_og = FRF_data.FRF_z;

f_res = length(f_vect)/f_vect(end);

f_max_interest = 8; % Hz
f_min_interest = 2;
n_f_max = ceil(f_max_interest*f_res); % number of elements of the first 8 Hz
n_f_min = ceil(f_min_interest*f_res); % number of elements of the first 2 Hz

f_vect = f_vect(n_f_min:n_f_max); % Frequency vector of interest
om_vect=f_vect*2*pi;
FRFs_og = FRFs_og(n_f_min:n_f_max,:)';  % FRFs of interest (2-8 Hz)
[m,n]=size(FRFs_og);

% look for FRFs in antimodal positions
%compute max value of each FRF to select only the stronger ones

%rough identfication of nodal positions
n_avg = 30;
[node,nodes_i] = findpeaks(mean((-1)*abs(FRFs_og(:,1:ceil(f_res*3.5))),1),"NPeaks",5,"MinPeakWidth",4);
n_modes = length(nodes_i)-1;
FRFs = zeros(n_modes,n_avg,n);
absavg = zeros(n_modes,n);
pks = zeros(n_modes,1);
pks_i = zeros(n_modes,1);
res_f = zeros(n_modes,1);
om1 = zeros(n_modes,1);
om2 = zeros(n_modes,1);
i1 = zeros(n_modes,1);
i2 = zeros(n_modes,1);
maxes_i_all = zeros(n_modes,m); 

for j = 1 : n_modes
    win = nodes_i(j):nodes_i(j+1);
    [~,rfi_rel] = max(mean(abs(FRFs_og(:,win)),1));   % rough resonance location inside this window
    rfi = win(rfi_rel);
    band = max(rfi-2,win(1)) : min(rfi+2,win(end));    % narrow band around that resonance
    maxes=zeros(1,m);
    for i=1:m
        maxes(i)=max(abs(imag(FRFs_og(i,band))))/mean(abs(imag(FRFs_og(i,:))));
    end
    %select first n_avg highest FRF around each fres
    [maxes,maxes_i]=sort(maxes,'descend');
    maxes_i_all(j,:)=maxes_i; %used later for nodes plots
    FRFi=maxes_i(1:n_avg);
    FRFs(j,:,:)=FRFs_og(FRFi,:);
    absavg(j,:)=mean(abs(FRFs(j,:,:)),2);

    % look for the single resonance peak of the mode restricted to its window
    [pks(j),pk_i_rel] = max(absavg(j,win));
    pks_i(j) = win(pk_i_rel);
    res_f(j) = f_vect(pks_i(j));

    %half power values
    i1(j)=find(absavg(j,1:pks_i(j))<absavg(j,pks_i(j))/sqrt(2),1,'last');
    i2(j)=find(absavg(j,pks_i(j):end)<absavg(j,pks_i(j))/sqrt(2),1,'first')+rfi;
    om1(j)=om_vect(i1(j));
    om2(j)=om_vect(i2(j));     
end
% h via half power method
om0 = res_f*2*pi;
h_vect=(om2.^2-om1.^2)./4./om0.^2;

%%
figure

for  j = 1 : n_modes
subplot(2,2,j)
semilogy(f_vect(1:ceil(3*f_res)), abs(squeeze(FRFs(j,:,1:ceil(3*f_res)))), 'cyan', LineWidth=.01)
hold on
semilogy(f_vect(1:ceil(3*f_res)),absavg(j,1:ceil(3*f_res)), 'magenta', LineWidth=2)
hold on
plot(f_vect(i1(j)),absavg(j,i1(j)),'o')
hold on
plot(f_vect(i2(j)),absavg(j,i2(j)),'o')
hold on
plot(res_f(j),pks(j),'bo',LineWidth=2)
grid on
xlim([2,5])


end
%%

%for each FRF the mode shape are scaled by a factor, fixed the input
% position in k, the mode shape  for FRF_kj relative to a i mode
% normalization of m and phi_input to 1
for i = 1:length(res_f)
    for j=1:m
        unnormed_mode(i,j)=-real(1i*2*om_vect(pks_i(i))^2*h_vect(i)*FRFs_og(j,pks_i(i)));
    end
end

load("FRF_H1.mat")
load("PuntiLaser_FS_FEM.mat")
load("modal_output.mat")
load("connectivity.mat")

% nodes in the Canopy from FE
[~,Loc1] = ismember(surface_canopy.Joint1,nodes.ID);
[~,Loc2] = ismember(surface_canopy.Joint2,nodes.ID);
[~,Loc3] = ismember(surface_canopy.Joint3,nodes.ID);
nodi123 = [Loc1,Loc2,Loc3];

% scale factor to display mode shapes
scalaFEM = 10;
scalaLASER = 70;

% plot FEM model + experimental measuring point seected as antinodal once per resonance frequency
%their location should give a first guess of the regions with highest displacement
for j = 1:n_modes
    figure(100+j);
    cc=ones(length(nodes.X),1);
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z],'CData',cc,'FaceColor','interp','EdgeColor','none')
    axis tight
    hold on;
    plot3(x(maxes_i_all(j,n_avg+1:end)),y(maxes_i_all(j,n_avg+1:end)),z(maxes_i_all(j,n_avg+1:end))+1,'ro','MarkerSize',4,'MarkerFaceColor','r')
    plot3(x(maxes_i_all(j,1:n_avg)),y(maxes_i_all(j,1:n_avg)),z(maxes_i_all(j,1:n_avg))+1,'go','MarkerSize',4,'MarkerFaceColor','g')
    title(sprintf('Selected antinodes for mode %d: f=%.3f Hz', j, res_f(j)))
end


for ii = 1:length(res_f)
    k=1;
    if ii==3
        k=-1;
    end

    %FEM
    mode = ii;
    mode_sel = modeshapes(modeshapes.No == mode,:); % This takes all rows of modeshapes elements of mode ii,
    % rows contain:  ux, uy, uz e ID
    [~,Locb] = ismember(mode_sel.ID,nodes.ID);
    modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]}; % scostamenti della modeshape ii in x y e z, dalla FEM

    % LASER - we only get z displacement
    fff = max(abs(modedef(:,3)));
    fff_exp=max(abs(k*unnormed_mode(ii,:)));
    F = scatteredInterpolant(x, y,k*unnormed_mode(ii,:)','natural'); % Funzione dei spostamenti laser registrati
    % calcolabile nei nodi della FEM
    F.ExtrapolationMethod = 'linear';


    z_LASER = F(nodes.X, nodes.Y);
    zPlot=z_LASER/fff_exp*scalaFEM; % ci sta

    figure
    cc0=z_LASER/fff_exp; % normalization
    ax1=subplot(1,2,1);
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y, nodes.Z+zPlot],...
        'CData',cc0,'FaceColor','interp','EdgeColor','none'); % pazzia
    title(sprintf('LASER -- Mode %i: f=%5.3f Hz', mode,res_f(ii)/sqrt(70)) )
    nmap = 10;
    map = jet(nmap);
    colormap(ax1,map)
    clim(ax1,[-1 1])
    cb0 = colorbar(ax1,'eastoutside');
    axis tight
    grid on
    view(2)
    %second o plot
    cc1=(modedef(:,3))/fff;
    ax2=subplot(1,2,2);
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z]+modedef/fff*scalaFEM,...
        'CData',cc1,'FaceColor','interp','EdgeColor','none');
    title(sprintf('FEM -- Mode %i: f=%5.3f Hz', mode,modpar.freq(mode)) )
    nmap1 = 10;
    map1 = jet(nmap1);
    colormap(ax2,map1)
    clim(ax2,[-1 1])
    cb = colorbar(ax2,'eastoutside');
    axis tight
    grid on
    view(2)

end

