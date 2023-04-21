function [] = test1()
n = 100;
x = linspace(-10,10,n);
y = linspace(-10,10,n);
z = zeros(10);
for i = 1:n 
    for j = 1:n
        z(i,j) = x(i)^2 + y(j)^2;
    end
end
contourf(x,y,z)
hold off
end

