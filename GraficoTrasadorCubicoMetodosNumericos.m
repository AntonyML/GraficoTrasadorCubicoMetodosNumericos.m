% Datos
x = [1, 1.5, 2];
f = [3, 2, 5];

% Calcular h
h = diff(x);

% Sistema tridiagonal para C
n = length(x) - 1;
A = zeros(n+1, n+1);
b = zeros(n+1, 1);

% Coeficientes de la matriz tridiagonal
A(1,1) = 1; % C0 = 0, natural spline
A(n+1, n+1) = 1; % Cn = 0, natural spline

for j = 2:n
    A(j, j-1) = h(j-1);
    A(j, j) = 2 * (h(j-1) + h(j));
    A(j, j+1) = h(j);
    b(j) = 3 * ((f(j+1) - f(j)) / h(j) - (f(j) - f(j-1)) / h(j-1));
end

% Resolver sistema
C = A \ b;

% Calcular b_j y d_j
b = zeros(n, 1);
d = zeros(n, 1);

for j = 1:n
    b(j) = (f(j+1) - f(j)) / h(j) - h(j) * (2*C(j) + C(j+1)) / 3;
    d(j) = (C(j+1) - C(j)) / (3 * h(j));
end

% Coeficientes a_j son los valores de f(x_j)
a = f;

% Evaluar S(x) en x = 1.25
x_eval = 1.25;
if x_eval < x(2)
    j = 1;
else
    j = 2;
end

S_eval = a(j) + b(j)*(x_eval - x(j)) + C(j)*(x_eval - x(j))^2 + d(j)*(x_eval - x(j))^3;
disp(['S(1.25) = ', num2str(S_eval)]);

% Graficar S(x)
xx = linspace(x(1), x(3), 100);
yy = zeros(size(xx));

for k = 1:length(xx)
    if xx(k) < x(2)
        j = 1;
    else
        j = 2;
    end
    yy(k) = a(j) + b(j)*(xx(k) - x(j)) + C(j)*(xx(k) - x(j))^2 + d(j)*(xx(k) - x(j))^3;
end

plot(x, f, 'ro', 'MarkerFaceColor', 'r'); % puntos de datos
hold on;
plot(xx, yy, 'b-'); % spline
xlabel('x');
ylabel('S(x)');
title('Trazador Cúbico');
grid on;
legend('Datos', 'Trazador Cúbico');
hold off;

