function D = distancias_pso(X, C)
%DISTANCIAS_PSO Calcula distancias euclidianas al cuadrado.
% Filas de X: pixeles; filas de C: centroides. D(i,k) es la distancia
% al cuadrado entre el pixel i y el centroide k.

    D = zeros(size(X, 1), size(C, 1));
    for k = 1:size(C, 1)
        diferencia = X - C(k, :);
        D(:, k) = sum(diferencia.^2, 2);
    end
end
