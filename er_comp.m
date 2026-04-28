function [err] = er_comp(Gnum,Gexp,x)
%Gnum è cell array di function handles, Gexp è un cell harray di valori della FRF di dimensione numero
%frequenze su cui ricostruisco FRF e numero punti dalle quali estraggo i
%dati [100,3]
[m,n]=size(Gnum);
err=[];
for i =1:m
    for j=1:n
        err=[err;(Gnum{i,j}(x)-Gexp{i,j})];
    end
end


end

