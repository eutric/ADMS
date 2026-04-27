function [outputArg1,outputArg2] = h_extim(Gexp,mode)
%take Gexp come valore sperimentale di FRF function handle
%mode come struct
om0=mode.OM;
omvect=linspace(0.8*om0,1.2*om0,100);
G_vect=Gexp(omvect);



end

