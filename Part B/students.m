clear
close all
clc

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
scalaLASER = 10;

% plot FEM model + experimental measuring point
figure(100);
cc=ones(length(nodes.X),1);
patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z],...
    'CData',cc,'FaceColor','interp','EdgeColor','none')
axis tight
hold on;
plot3(x,y,z+1,'ro','MarkerSize',4,'MarkerFaceColor','r')


%%

for ii = 1:4

    %FEM
    mode = ii;
    mode_sel = modeshapes(modeshapes.No == mode,:);
    [~,Locb] = ismember(mode_sel.ID,nodes.ID);
    modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]};

    fff = max(abs(modedef(:,3)));

    figure(ii);
    axis equal
    grid on
    hold on
    view(2)

    cc=(modedef(:,3))/fff;
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z]+modedef/fff*scalaFEM,...
        'CData',cc,'FaceColor','interp','EdgeColor','none')
    title(sprintf('FEM -- Mode %i: f=%5.3f Hz', mode,modpar.freq(mode)) )
    nmap = 10;
    map     = jet(nmap);
    colormap(gca,map)
    clim([-1 1])
    cb = colorbar('eastoutside');
    axis tight

end