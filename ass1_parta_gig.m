%countilever beam parameters
clear all
close all
clc
L=1200;  % mm
h=8;
b=40;
rho=2700e-9; %kg/mm^3
E=68000; %MPa
m=rho*h*b; % kg/mm
J=h^3*b/12; % mm^4;
f_vect=linspace(0,10,10000);
om_vect=2*pi.*f_vect;
g=@(om)(m*om.^2/E/J).^(1/4);

%build the H matrix
H=@(g)[1 0 1 0;
    0 g 0 g;
    -E*J.*g.^2.*cos(g*L) -E*J.*g.^2.*sin(g*L) E*J.*g.^2.*cosh(g*L) E*J.*g.^2.*sinh(g*L);
    E*J.*g.^3.*sin(g*L) -E*J*g.^3.*cos(g*L) E*J*g.^3.*sinh(g*L) E*J*g.^3.*cosh(g*L) ];
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
%natural frequencies
%computing mode shapes
figure
xx=linspace(0,1200,1000);
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


