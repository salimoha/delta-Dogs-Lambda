function [ y] = test_fun( x )

y = sin(30.*(x-0.9).^4).*cos(2.*(x-0.9)) + (x-0.9)./2;


end

