function y = rastriginn(x)

global n

A=2;
y = A * n * ones(1,size(x,2));
for ii = 1 : 1 : n
    y = y + (x(ii,:) .^ 2 - A * cos(2 * pi * x(ii,:)));
end

end
