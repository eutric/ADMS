function [h] = h_extim(Gexp,om0)
%take Gexp come valore sperimentale di FRF function handle
%mode come struct

omvect=linspace(0.9*om0,1.1*om0,100);
G_vect=Gexp(omvect);
m=max(abs(G_vect));
i_om1=find(abs(G_vect)>m/sqrt(2),1);
i_om2=find(abs(G_vect)>m/sqrt(2),1,"last");
om1=omvect(i_om1);
om2=omvect(i_om2);
h=(om2.^2-om1.^2)/4/om0/om0;

end

