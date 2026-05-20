function c = mid_node(a, b, p)
    if nargin == 2
        ac = 1/2*(b-a);
    else
        ac = p*(b-a);
    end
    c = a + ac;
end