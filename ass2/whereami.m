function el = whereami(P, xy)
% This now works only for elements that don't cross each other x
% coordinates, they have to go from  left to right (o viceversa)
el = [];
for ii=1:length(xy(:,1))-1
    if P(1) >= xy(ii,1) && P(1) < xy(ii+1,1)
        if P(2) >= xy(ii,2) 
            if abs(xy(ii,2)-xy(ii+1,2)) < 1e-6
                el = ii;
            elseif (P(2) < xy(ii+1,2))
                el = ii;
            end
        end
    end
end
if isempty(el)
    el = ii;
end
end