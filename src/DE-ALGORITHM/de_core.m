
function resultado = de_core(Xm, parametros)
% DE_CORE Ejecuta Evolucion Diferencial para optimizar centroides.
%
% Variante implementada: DE/rand/1/bin.
%
% Entradas:
%   Xm         - Muestra normalizada de pixeles.
%   parametros - Estructura con los parametros de DE.
%
% Campos de parametros:
%   K             - Numero de centroides.
%   NP            - Tamano de la poblacion.
%   F             - Factor de mutacion.
%   CR            - Probabilidad de cruce.
%   maxIter       - Numero maximo de generaciones.
%   penalizacion  - Penalizacion por grupos vacios.
%
% Salida:
%   resultado.mejoresCentros
%   resultado.mejorCosteGlobal
%   resultado.historial
%   resultado.tiempo
%   resultado.mejorIndividuo
%   resultado.costeInicial
%
% Cada individuo representa K centroides.
% Para K=4 y tres indices espectrales existen 12 variables.
%
% La funcion objetivo es la misma que utilizan PSO y SA.

%% 1. Obtener parametros

K = parametros.K;
NP = parametros.NP;
F = parametros.F;
CR = parametros.CR;
maxIter = parametros.maxIter;
penalizacion = parametros.penalizacion;

%% 2. Comprobar parametros

assert(isnumeric(Xm) && ismatrix(Xm) && ~isempty(Xm), ...
    'Xm debe ser una matriz numerica no vacia.');

assert(all(isfinite(Xm(:))), ...
    'Xm contiene valores no finitos.');

assert(K >= 1 && K == floor(K), ...
    'K debe ser un entero positivo.');

assert(NP >= 4 && NP == floor(NP), ...
    'NP debe ser un entero mayor o igual que 4.');

assert(maxIter >= 1 && maxIter == floor(maxIter), ...
    'maxIter debe ser un entero positivo.');

assert(F > 0 && F <= 2, ...
    'F debe estar en el intervalo (0, 2].');

assert(CR >= 0 && CR <= 1, ...
    'CR debe estar en el intervalo [0, 1].');

assert(penalizacion >= 0, ...
    'La penalizacion no puede ser negativa.');

nDimensiones = size(Xm, 2);

assert(size(Xm, 1) >= K, ...
    'No hay suficientes pixeles para inicializar los centroides.');

%% 3. Definir el espacio de busqueda

% Cada individuo contiene K centroides.
nVariables = K * nDimensiones;

% Limites obtenidos de la muestra de optimizacion.
minX = min(Xm, [], 1);
maxX = max(Xm, [], 1);

rango = maxX - minX;


limiteInferior = repmat(minX, 1, K);
limiteSuperior = repmat(maxX, 1, K);

%% 4. Inicializar la poblacion

% Cada fila representa un individuo.
% Cada individuo contiene K centroides.

poblacion = zeros(NP, nVariables);

for i = 1:NP

    % Inicializacion aleatoria uniforme
    % dentro de los limites de cada indice.

    centros = zeros(K, nDimensiones);

    for k = 1:K

        centros(k, :) = minX + rand(1, nDimensiones) ...
            .* rango;

    end

    poblacion(i, :) = reshape(centros.', 1, []);

end

%% 5. Evaluar la poblacion inicial

costes = zeros(NP, 1);

for i = 1:NP

    centros = reshape( ...
        poblacion(i, :), nDimensiones, K).';

    costes(i) = calcular_coste( ...
        Xm, centros, penalizacion);

end

% Mejor individuo inicial.
[mejorCosteGlobal, mejorIndice] = min(costes);

mejorIndividuo = poblacion(mejorIndice, :);

costeInicial = mejorCosteGlobal;

% Historial del mejor coste por generacion.
historial = zeros(maxIter, 1);

%% 6. Ejecutar Evolucion Diferencial

tic;

for iter = 1:maxIter

    % Crear la nueva poblacion.
    nuevaPoblacion = poblacion;

    nuevosCostes = costes;

    for i = 1:NP

        %% 6.1 Seleccionar tres individuos diferentes

        % Los individuos seleccionados deben ser:
        % diferentes entre si y distintos del individuo i.

        candidatos = setdiff(1:NP, i);

        indices = randperm(NP - 1, 3);

        r1 = candidatos(indices(1));
        r2 = candidatos(indices(2));
        r3 = candidatos(indices(3));

        %% 6.2 Mutacion diferencial

        % Variante DE/rand/1:
        %
        % v = x_r1 + F * (x_r2 - x_r3)

        vectorMutante = poblacion(r1, :) ...
            + F * (poblacion(r2, :) - poblacion(r3, :));

        %% 6.3 Control de limites

        % Mantener las variables dentro del espacio
        % de busqueda definido por los datos.

        vectorMutante = max( ...
            min(vectorMutante, limiteSuperior), ...
            limiteInferior);

        %% 6.4 Cruce binomial

        vectorPrueba = poblacion(i, :);

        % Garantizar al menos una componente procedente
        % del vector mutante.

        jRand = randi(nVariables);

        for j = 1:nVariables

            if rand <= CR || j == jRand

                vectorPrueba(j) = vectorMutante(j);

            end

        end

        %% 6.5 Evaluar el individuo de prueba

        centrosPrueba = reshape( ...
            vectorPrueba, nDimensiones, K).';

        costePrueba = calcular_coste( ...
            Xm, centrosPrueba, penalizacion);

        %% 6.6 Seleccion

        % Se conserva el individuo con menor coste.

        if costePrueba <= costes(i)

            nuevaPoblacion(i, :) = vectorPrueba;

            nuevosCostes(i) = costePrueba;

        end

    end

    %% 6.7 Actualizar la poblacion

    poblacion = nuevaPoblacion;

    costes = nuevosCostes;

    %% 6.8 Actualizar la mejor solucion global

    [mejorCosteGeneracion, mejorIndice] = min(costes);

    if mejorCosteGeneracion < mejorCosteGlobal

        mejorCosteGlobal = mejorCosteGeneracion;

        mejorIndividuo = poblacion(mejorIndice, :);

    end

    %% 6.9 Registrar convergencia

    historial(iter) = mejorCosteGlobal;

    if mod(iter, 20) == 0 || iter == 1

        fprintf(['Generacion %d/%d | ', ...
            'Mejor coste: %.6f\n'], ...
            iter, maxIter, mejorCosteGlobal);

    end

end

tiempo = toc;

%% 7. Recuperar los mejores centroides

mejoresCentros = reshape( ...
    mejorIndividuo, nDimensiones, K).';

%% 8. Preparar los resultados

resultado.mejoresCentros = mejoresCentros;

resultado.mejorCosteGlobal = mejorCosteGlobal;

resultado.historial = historial;

resultado.tiempo = tiempo;

resultado.mejorIndividuo = mejorIndividuo;

resultado.costeInicial = costeInicial;

resultado.poblacionFinal = poblacion;

resultado.costesFinales = costes;

end

%% FUNCION AUXILIAR: CALCULAR COSTE

function coste = calcular_coste(X, centros, penalizacion)

    % Calcular distancias al cuadrado.
    D = distancias_de(X, centros);

    % Obtener distancia minima de cada pixel.
    [distanciasMinimas, etiquetas] = min(D, [], 2);

    % Coste base: distancia media al cuadrado.
    costeBase = mean(distanciasMinimas);

    % Contar grupos vacios.
    gruposVacios = size(centros, 1) ...
        - numel(unique(etiquetas));

    % Funcion objetivo penalizada.
    coste = costeBase + penalizacion * gruposVacios;

end
