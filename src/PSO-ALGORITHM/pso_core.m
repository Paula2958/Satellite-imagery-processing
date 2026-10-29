function resultado = pso_core(Xm, parametros)
% PSO_CORE Ejecuta Particle Swarm Optimization para agrupar datos.
%
% Entradas:
%   Xm         - Muestra normalizada de pixeles e indices espectrales.
%   parametros - Estructura con los parametros de PSO.
%
% Salida:
%   resultado  - Mejor solucion, coste e historial de convergencia.

%% 1. Obtener parametros
K = parametros.K;
nParticulas = parametros.nParticulas;
maxIter = parametros.maxIter;
w = parametros.w;
c1 = parametros.c1;
c2 = parametros.c2;
penalizacion = parametros.penalizacion;

nDimensiones = size(Xm, 2);
nVariables = K * nDimensiones;

%% 2. Definir limites de los centroides
minX = min(Xm, [], 1);
maxX = max(Xm, [], 1);

rango = maxX - minX;
rango(rango == 0) = 1;

limiteInferior = repmat(minX, 1, K);
limiteSuperior = repmat(maxX, 1, K);

velocidadMaxima = 0.2 * repmat(rango, 1, K);

%% 3. Inicializar las particulas
% Cada particula representa K centroides.
posiciones = zeros(nParticulas, nVariables);
velocidades = zeros(nParticulas, nVariables);

for p = 1:nParticulas
    centros = zeros(K, nDimensiones);

    % El primer centro se elige aleatoriamente.
    primerIndice = randi(size(Xm, 1));
    centros(1, :) = Xm(primerIndice, :);

    % Los siguientes centros se eligen buscando puntos
    % alejados de los centros seleccionados.
    for k = 2:K
        D = distancias_pso(Xm, centros(1:k-1, :));
        distanciaMinima = min(D, [], 2);

        [~, indice] = max(distanciaMinima);
        centros(k, :) = Xm(indice, :);
    end

    posiciones(p, :) = reshape(centros.', 1, []);
end

%% 4. Evaluar las posiciones iniciales
mejoresPosiciones = posiciones;
mejoresCostes = inf(nParticulas, 1);

for p = 1:nParticulas
    centros = reshape(posiciones(p, :), nDimensiones, K).';
    mejoresCostes(p) = calcular_coste(Xm, centros, penalizacion);
end

[mejorCosteGlobal, mejorIndice] = min(mejoresCostes);
mejorPosicionGlobal = mejoresPosiciones(mejorIndice, :);

historial = zeros(maxIter, 1);

%% 5. Bucle principal de PSO
tic;

for iter = 1:maxIter
    for p = 1:nParticulas

        r1 = rand(1, nVariables);
        r2 = rand(1, nVariables);

        % Actualizar la velocidad de cada particula.
        velocidades(p, :) = w * velocidades(p, :) ...
            + c1 * r1 .* ...
            (mejoresPosiciones(p, :) - posiciones(p, :)) ...
            + c2 * r2 .* ...
            (mejorPosicionGlobal - posiciones(p, :));

        % Limitar la velocidad.
        velocidades(p, :) = max( ...
            min(velocidades(p, :), velocidadMaxima), ...
            -velocidadMaxima);

        % Actualizar la posicion de la particula.
        posiciones(p, :) = posiciones(p, :) + velocidades(p, :);

        % Mantener los centroides dentro de los limites.
        posiciones(p, :) = max( ...
            min(posiciones(p, :), limiteSuperior), ...
            limiteInferior);

        centros = reshape(posiciones(p, :), nDimensiones, K).';

        costeActual = calcular_coste(Xm, centros, penalizacion);

        % Actualizar la mejor solucion individual.
        if costeActual < mejoresCostes(p)
            mejoresCostes(p) = costeActual;
            mejoresPosiciones(p, :) = posiciones(p, :);
        end

        % Actualizar la mejor solucion global.
        if mejoresCostes(p) < mejorCosteGlobal
            mejorCosteGlobal = mejoresCostes(p);
            mejorPosicionGlobal = mejoresPosiciones(p, :);
        end
    end

    historial(iter) = mejorCosteGlobal;

    if mod(iter, 20) == 0 || iter == 1
        fprintf('Iteracion %d/%d | mejor coste: %.6f\n', ...
            iter, maxIter, mejorCosteGlobal);
    end
end

tiempo = toc;

%% 6. Preparar los resultados
mejoresCentros = reshape( ...
    mejorPosicionGlobal, nDimensiones, K).';

resultado.mejoresCentros = mejoresCentros;
resultado.mejorCosteGlobal = mejorCosteGlobal;
resultado.historial = historial;
resultado.tiempo = tiempo;

end

%% Funcion auxiliar: calcular el coste
function coste = calcular_coste(X, centros, penalizacion)

D = distancias_pso(X, centros);
distanciasMinimas = min(D, [], 2);

costeBase = mean(distanciasMinimas);

[~, etiquetas] = min(D, [], 2);
gruposVacios = size(centros, 1) - numel(unique(etiquetas));

coste = costeBase + penalizacion * gruposVacios;

end
