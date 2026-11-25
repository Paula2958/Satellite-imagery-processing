
function D = distancias_de(X, C)
%DISTANCIAS_DE Calcula distancias euclidianas al cuadrado.
%
% Entradas:
%   X - Matriz de datos. Cada fila representa un pixel.
%   C - Matriz de centroides. Cada fila representa un centroide.
%
% Salida:
%   D - Matriz de distancias al cuadrado.
%       D(i,k) representa la distancia entre X(i,:) y C(k,:).
%
% Esta implementacion mantiene la misma metrica utilizada
% por los algoritmos PSO y SA.

D = zeros(size(X, 1), size(C, 1));

for k = 1:size(C, 1)

    diferencia = X - C(k, :);

    D(:, k) = sum(diferencia.^2, 2);

end

end
