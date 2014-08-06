function f=fact(order);
% FACT: Factorial function
f = 1;
if order > 0;
        for j=1:order;
        f = f .* j;
        end
        end

