function [err] = er_comp(Gnum,Gexp,x)
%Gnum è cell array di function handles, Gexp è un cell harray di valori della FRF di dimensione numero
%frequenze su cui ricostruisco FRF e numero punti dalle quali estraggo i
%dati [100,3]
err=0;
[m,n]=size(Gnum);
err=0;
for i =1:m
    for j=1:n
        err=err+real(Gnum{i,j}(x)-Gexp{i,j}).^2+imag(Gnum{i,j}(x)-Gexp{i,j}).^2;
    end
end


end

