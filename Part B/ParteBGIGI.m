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
FRFs_og = FRFs_og(n_f_min:n_f_max,:)';  % FRFs of interest (0-8 Hz)
[m,n]=size(FRFs_og);

% look for FRFs in antimodal positions
%compute max value of each FRF to select only the stronger ones
max_v=[];
for i=1:m
    max_v(i)=max(FRFs_og(i,:));
end
[maxes,maxes_i]=sort(max_v,'descend');
%select first 50 FRF
FRFi=maxes_i(1:100);
FRFs=FRFs_og(FRFi,:);
absavg=mean(abs(FRFs),1);
phavg=mean(angle(FRFs),1);
figure
subplot(2,1,1)
semilogy(f_vect, abs(FRFs), 'cyan', LineWidth=.01)
grid on
xlim([2,8])
hold on
semilogy(f_vect, absavg, 'magenta', LineWidth=2)

subplot(2,1,2)
plot(f_vect, angle(FRFs), 'cyan',LineWidth=.01)
xlim([2,8])
grid on


%% look fer peaks to identify the resonant frequencies

[pks,pks_i]=findpeaks(absavg,"MinPeakWidth",5,"Threshold",2.05*1e-5);%using this criterias only the relevant peaks are found
subplot(2,1,1)
plot(f_vect(pks_i),pks,'bo',LineWidth=2)
res_f=f_vect(pks_i);
for i = 1:length(res_f)
    rf=res_f(i);%resonant frequency
    rfi=pks_i(i);%index of resonsnt freq in f_vect
    i1=find(absavg(1:rfi)<absavg(rfi)/sqrt(2),1,'last');
    subplot(2,1,1)
    plot(f_vect(i1),absavg(i1),'o')
    i2=find(absavg(rfi:end)<absavg(rfi)/sqrt(2),1,'first')+rfi;
    subplot(2,1,1)
    plot(f_vect(i2),absavg(i2),'o')
    % compute f1 and f2
    om1(i)=om_vect(i1);
    om2(i)=om_vect(i2);
    om0(i)=om_vect(rfi);
end
h_vect=(om2.^2-om1.^2)./4./om0.^2;
%for each FRF the mode shape are scaled of a factor, fixed the input
%position in k, the mode shape  for FRF_kj relative to a i mode
for i =1:length(res_f)
    for j=1:m
        unnormed_mode(i,j)=-imag(2*h_vect(i).*(om_vect(pks_i(i))).^2*FRFs_og(j,pks_i(i)));
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

% plot FEM model + experimental measuring point
figure(100);
cc=ones(length(nodes.X),1);
patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z],...
    'CData',cc,'FaceColor','interp','EdgeColor','none')
axis tight
hold on;
plot3(x,y,z+1,'ro','MarkerSize',4,'MarkerFaceColor','r')


for ii = 1:4

    %FEM
    mode = ii;
    mode_sel = modeshapes(modeshapes.No == mode,:); % This takes all rows of modeshapes elements of mode ii,
    % rows contain:  ux, uy, uz e ID
    [~,Locb] = ismember(mode_sel.ID,nodes.ID);
    modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]}; % scostamenti della modeshape ii in x y e z, dalla FEM

    % LASER - we only get z displacement
    fff = max(abs(modedef(:,3)));

    F = scatteredInterpolant(x, y, scalaLASER/fff*unnormed_mode(ii,:)','natural'); % Funzione dei spostamenti laser registrati
    Z = scatteredInterpolant(x,y,z,'natural');
    % calcolabile nei nodi della FEM

    z_LASER = F(nodes.X, nodes.Y);
    z_nodes=Z(nodes.X, nodes.Y);

    figure 
    axis equal
    grid on
    hold on
    view(2)

    cc=z_LASER/fff; % normalization
    figure(ii)
    subplot(1,2,1)
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y, z_LASER],...
        'CData',cc,'FaceColor','interp','EdgeColor','none');
    title(sprintf('FEM -- Mode %i: f=%5.3f Hz', mode,modpar.freq(mode)) )
    nmap = 10;
    map     = jet(nmap);
    colormap(gca,map)
    clim([-1 1])
    cb = colorbar('eastoutside');
    axis tight
    subplot(1,2,2)
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z-z_nodes]+modedef/fff*scalaFEM,...
        'CData',cc,'FaceColor','interp','EdgeColor','none');
    title(sprintf('FEM -- Mode %i: f=%5.3f Hz', mode,modpar.freq(mode)) )
    nmap1 = 10;
    map1     = jet(nmap1);
    colormap(gcf,map1)
    clim([-1 1])
    cb = colorbar('eastoutside');
    axis tight

end




