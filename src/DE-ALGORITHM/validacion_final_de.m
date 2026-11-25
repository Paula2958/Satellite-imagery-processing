
%% VALIDACION_FINAL_DE
% Validacion interna de la mejor solucion de Evolucion Diferencial.
%
% Comprueba:
%   - Coherencia del mapa de segmentacion.
%   - Reproducibilidad de la funcion objetivo.
%   - Distribucion y compactacion de los grupos.
%   - Medias espectrales reales.
%   - Separacion entre centroides.
%   - Indice de Davies-Bouldin.
%   - Historial de convergencia.
%
% No vuelve a ejecutar la optimizacion.
% No requiere toolboxes adicionales.

clearvars;
clc;
close all;

%% 1. Definir rutas

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));

rawDataDir = fullfile(projectRoot, 'data', 'raw');
resultsDir = fullfile(projectRoot, 'results');
figuresDir = fullfile(resultsDir, 'figures', 'DE');

addpath(scriptDir);

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

fprintf('============================================\n');
fprintf(' VALIDACION FINAL DE EVOLUCION DIFERENCIAL\n');
fprintf('============================================\n');

%% 2. Cargar la mejor solucion

archivo = fullfile(resultsDir, 'de_mejor_resultado.mat');

assert(exist(archivo, 'file') == 2, ...
    'Ejecuta primero de_mejor_resultado.m.');

datos = load(archivo);

parametros = datos.parametros;
K = parametros.K;

mejoresCentros = datos.mejoresCentros;
centrosOriginales = datos.centrosOriginales;

mapaClases = datos.mapaClases;
etiquetas = datos.etiquetas;
historial = datos.mejorHistorial;

fprintf('Semilla seleccionada: %d\n', datos.mejorSemilla);
fprintf('Coste de optimizacion: %.6f\n', datos.mejorCoste);

%% 3. Cargar bandas Sentinel-2

B02 = cargar_banda_de(rawDataDir, '*_B02_(Raw).tiff');
B03 = cargar_banda_de(rawDataDir, '*_B03_(Raw).tiff');
B04 = cargar_banda_de(rawDataDir, '*_B04_(Raw).tiff');
B08 = cargar_banda_de(rawDataDir, '*_B08_(Raw).tiff');
B11 = cargar_banda_de(rawDataDir, '*_B11_(Raw).tiff');

assert(isequal(size(B02), size(B03), size(B04), ...
    size(B08), size(B11)), ...
    'Las bandas no tienen las mismas dimensiones.');

%% 4. Reconstruir mascara de pixeles validos

mascara = isfinite(B02) & isfinite(B03) & ...
          isfinite(B04) & isfinite(B08) & ...
          isfinite(B11) & ...
          B02 > 0 & B03 > 0 & B04 > 0 & ...
          B08 > 0 & B11 > 0;

%% 5. Reconstruir indices espectrales

denNDVI = B08 + B04;
denNDWI = B03 + B08;
denNDBI = B11 + B08;

NDVI = nan(size(B02));
NDWI = nan(size(B02));
NDBI = nan(size(B02));

validosNDVI = mascara & denNDVI ~= 0;
validosNDWI = mascara & denNDWI ~= 0;
validosNDBI = mascara & denNDBI ~= 0;

NDVI(validosNDVI) = ...
    (B08(validosNDVI) - B04(validosNDVI)) ./ ...
    denNDVI(validosNDVI);

NDWI(validosNDWI) = ...
    (B03(validosNDWI) - B08(validosNDWI)) ./ ...
    denNDWI(validosNDWI);

NDBI(validosNDBI) = ...
    (B11(validosNDBI) - B08(validosNDBI)) ./ ...
    denNDBI(validosNDBI);

%% 6. Reconstruir matriz de datos

pixelesMascara = find(mascara);

Xtotal = [NDVI(mascara), ...
          NDWI(mascara), ...
          NDBI(mascara)];

filasFinitas = all(isfinite(Xtotal), 2);

X = Xtotal(filasFinitas, :);

pixelesValidos = pixelesMascara(filasFinitas);

nPixeles = size(X, 1);
nPixelesTotales = numel(B02);
nPixelesExcluidos = nPixelesTotales - nPixeles;

assert(nPixeles > 0, ...
    'No hay pixeles validos.');

% Utilizar la normalizacion guardada en la mejor solucion.

Xnormalizado = (X - datos.media) ./ datos.desv;

tolerancia = 1e-10;

%% 7. Validacion estructural

fprintf('\n--- VALIDACION ESTRUCTURAL ---\n');

assert(isequal(size(mapaClases), size(B02)), ...
    'El mapa no coincide con las dimensiones de la imagen.');

assert(numel(etiquetas) == nPixeles, ...
    'El numero de etiquetas no coincide con los pixeles validos.');

assert(isequal(pixelesValidos(:), datos.pixelesValidos(:)), ...
    'Los indices de pixeles validos no coinciden.');

assert(all(isfinite(etiquetas(:))), ...
    'Existen etiquetas no finitas.');

assert(all(etiquetas(:) >= 1 & etiquetas(:) <= K), ...
    'Existen etiquetas fuera del intervalo permitido.');

assert(all(etiquetas(:) == floor(etiquetas(:))), ...
    'Existen etiquetas no enteras.');

assert(isequal(size(mejoresCentros), [K, size(X, 2)]), ...
    'Las dimensiones de los centroides son incorrectas.');

assert(all(isfinite(mejoresCentros(:))), ...
    'Los centroides contienen valores no finitos.');

% Verificar posicion exacta de los pixeles clasificados.

mapaEsperado = zeros(size(B02), 'uint8');
mapaEsperado(pixelesValidos) = uint8(etiquetas);

assert(isequal(mapaClases, mapaEsperado), ...
    'El mapa no coincide con las etiquetas guardadas.');

fprintf('Dimensiones correctas: %d x %d\n', ...
    size(mapaClases, 1), size(mapaClases, 2));

fprintf('Pixeles totales: %d\n', nPixelesTotales);
fprintf('Pixeles clasificados: %d\n', nPixeles);
fprintf('Pixeles excluidos: %d\n', nPixelesExcluidos);

fprintf('Porcentaje clasificado: %.2f%%\n', ...
    100 * nPixeles / nPixelesTotales);

fprintf('Validacion estructural: CORRECTA\n');

%% 8. Validacion de la funcion objetivo

fprintf('\n--- VALIDACION NUMERICA ---\n');

D = distancias_de(Xnormalizado, mejoresCentros);

[distanciasMinimas, etiquetasVerificadas] = min(D, [], 2);

assert(isequal(etiquetas(:), etiquetasVerificadas(:)), ...
    'Las etiquetas no coinciden con el centroide mas cercano.');

conteo = accumarray(etiquetasVerificadas, 1, [K, 1]);

gruposVacios = sum(conteo == 0);

costeBase = mean(distanciasMinimas);

costeValidado = costeBase + ...
    parametros.penalizacion * gruposVacios;

diferenciaCoste = abs(costeValidado - datos.costeTotal);

assert(diferenciaCoste < tolerancia, ...
    'El coste recalculado no coincide con el coste guardado.');

assert(isequal(conteo(:), datos.conteoTotal(:)), ...
    'El conteo de grupos no coincide con el guardado.');

assert(gruposVacios == 0, ...
    'La solucion contiene grupos vacios.');

fprintf('Coste guardado: %.6f\n', datos.costeTotal);
fprintf('Coste recalculado: %.6f\n', costeValidado);
fprintf('Diferencia de coste: %.12f\n', diferenciaCoste);
fprintf('Grupos vacios: %d\n', gruposVacios);

fprintf('Validacion numerica: CORRECTA\n');

%% 9. Estadisticas espectrales por grupo

fprintf('\n--- VALIDACION ESPECTRAL ---\n');

mediasEspectrales = zeros(K, 3);
desvEspectrales = zeros(K, 3);

mediasNormalizadas = zeros(K, 3);

compactacion = zeros(K, 1);
dispersionMedia = zeros(K, 1);

desplazamientoCentroide = zeros(K, 1);

porcentajes = 100 * conteo / nPixeles;

for k = 1:K

    pertenece = etiquetasVerificadas == k;

    Xgrupo = X(pertenece, :);
    XgrupoNormalizado = Xnormalizado(pertenece, :);

    % Estadisticas reales de los pixeles asignados.

    mediasEspectrales(k, :) = mean(Xgrupo, 1);

    desvEspectrales(k, :) = std(Xgrupo, 0, 1);

    mediasNormalizadas(k, :) = mean(XgrupoNormalizado, 1);

    % Compactacion respecto al centroide optimizado.
    % Distancia euclidiana media al cuadrado.

    diferencias = XgrupoNormalizado - mejoresCentros(k, :);

    compactacion(k) = mean(sum(diferencias.^2, 2));

    % Dispersion respecto a la media real del grupo.
    % Distancia euclidiana media para Davies-Bouldin.

    diferenciasMedia = ...
        XgrupoNormalizado - mediasNormalizadas(k, :);

    dispersionMedia(k) = ...
        mean(sqrt(sum(diferenciasMedia.^2, 2)));

    % Diferencia entre centroide optimizado y media real.
    % Se calcula en el espacio normalizado.

    desplazamientoCentroide(k) = norm( ...
        mejoresCentros(k, :) - mediasNormalizadas(k, :));

end

tablaEspectral = table( ...
    (1:K)', ...
    conteo, ...
    porcentajes, ...
    mediasEspectrales(:, 1), ...
    mediasEspectrales(:, 2), ...
    mediasEspectrales(:, 3), ...
    desvEspectrales(:, 1), ...
    desvEspectrales(:, 2), ...
    desvEspectrales(:, 3), ...
    compactacion, ...
    desplazamientoCentroide, ...
    'VariableNames', { ...
        'Grupo', ...
        'Pixeles', ...
        'Porcentaje', ...
        'NDVI_Medio', ...
        'NDWI_Medio', ...
        'NDBI_Medio', ...
        'NDVI_Desv', ...
        'NDWI_Desv', ...
        'NDBI_Desv', ...
        'Compactacion', ...
        'DesplazamientoCentroide'});

disp(tablaEspectral);

% Comprobar que la compactacion ponderada reproduce el coste base.

compactacionGlobal = sum(conteo .* compactacion) / nPixeles;

assert(abs(compactacionGlobal - costeBase) < tolerancia, ...
    'La compactacion global no coincide con el coste base.');

fprintf('Compactacion global: %.6f\n', compactacionGlobal);

%% 10. Separacion entre centroides optimizados

fprintf('\n--- SEPARACION ENTRE CENTROIDES ---\n');

matrizSeparacion = zeros(K, K);

for i = 1:K

    for j = 1:K

        matrizSeparacion(i, j) = norm( ...
            mejoresCentros(i, :) - mejoresCentros(j, :));

    end

end

% Considerar cada pareja de centroides una sola vez.

mascaraPares = triu(true(K), 1);

separaciones = matrizSeparacion(mascaraPares);

separacionMinima = min(separaciones);
separacionMedia = mean(separaciones);
separacionMaxima = max(separaciones);

fprintf('Matriz de distancias entre centroides:\n');

disp(matrizSeparacion);

fprintf('Separacion minima: %.6f\n', separacionMinima);
fprintf('Separacion media: %.6f\n', separacionMedia);
fprintf('Separacion maxima: %.6f\n', separacionMaxima);

%% 11. Indice de Davies-Bouldin

fprintf('\n--- INDICE DE DAVIES-BOULDIN ---\n');

% Se utilizan las medias reales de cada grupo para calcular
% el indice Davies-Bouldin convencional.
%
% Rij = (Si + Sj) / Mij
%
% Si: distancia media de los pixeles a la media del grupo i.
% Mij: distancia entre las medias de los grupos i y j.

R = -inf(K, K);

for i = 1:K

    for j = 1:K

        if i ~= j

            distancia = norm( ...
                mediasNormalizadas(i, :) - ...
                mediasNormalizadas(j, :));

            if distancia <= eps

                R(i, j) = Inf;

            else

                R(i, j) = ...
                    (dispersionMedia(i) + ...
                     dispersionMedia(j)) / distancia;

            end

        end

    end

end

indiceDB = mean(max(R, [], 2));

fprintf('Indice Davies-Bouldin: %.6f\n', indiceDB);
fprintf('Un valor menor indica mejor calidad interna.\n');

if ~isfinite(indiceDB)

    warning('El indice Davies-Bouldin no es finito.');

end

%% 12. Validar historial de convergencia

fprintf('\n--- VALIDACION DE CONVERGENCIA ---\n');

assert(numel(historial) == parametros.maxIter, ...
    'La longitud del historial no es correcta.');

assert(all(isfinite(historial(:))), ...
    'El historial contiene valores no finitos.');

convergenciaMonotona = ...
    all(diff(historial(:)) <= tolerancia);

assert(convergenciaMonotona, ...
    'El historial presenta incrementos en el mejor coste.');

assert(abs(historial(end) - datos.mejorCoste) < tolerancia, ...
    'El ultimo coste no coincide con el mejor coste guardado.');

costePrimeraGeneracion = historial(1);
costeFinal = historial(end);

reduccionCoste = 100 * ...
    (costePrimeraGeneracion - costeFinal) / ...
    costePrimeraGeneracion;

fprintf('Coste primera generacion: %.6f\n', ...
    costePrimeraGeneracion);

fprintf('Coste ultima generacion: %.6f\n', costeFinal);

fprintf('Reduccion del coste: %.2f%%\n', reduccionCoste);

fprintf('Historial no creciente: CORRECTO\n');

%% 13. Guardar tablas de validacion

writetable(tablaEspectral, ...
    fullfile(resultsDir, 'de_validacion_clases.csv'));

tablaValidacion = table( ...
    datos.mejorSemilla, ...
    datos.mejorCoste, ...
    costeValidado, ...
    datos.mejorTiempo, ...
    nPixeles, ...
    nPixelesExcluidos, ...
    gruposVacios, ...
    compactacionGlobal, ...
    separacionMinima, ...
    separacionMedia, ...
    separacionMaxima, ...
    indiceDB, ...
    reduccionCoste, ...
    convergenciaMonotona, ...
    'VariableNames', { ...
        'Semilla', ...
        'CosteMuestra', ...
        'CosteTotal', ...
        'Tiempo', ...
        'PixelesValidos', ...
        'PixelesExcluidos', ...
        'GruposVacios', ...
        'CompactacionGlobal', ...
        'SeparacionMinima', ...
        'SeparacionMedia', ...
        'SeparacionMaxima', ...
        'DaviesBouldin', ...
        'ReduccionCoste', ...
        'ConvergenciaMonotona'});

writetable(tablaValidacion, ...
    fullfile(resultsDir, 'de_validacion_resumen.csv'));

%% 14. Guardar resultados completos

save(fullfile(resultsDir, 'de_validacion.mat'), ...
    'tablaEspectral', ...
    'tablaValidacion', ...
    'mediasEspectrales', ...
    'desvEspectrales', ...
    'mediasNormalizadas', ...
    'compactacion', ...
    'dispersionMedia', ...
    'compactacionGlobal', ...
    'desplazamientoCentroide', ...
    'matrizSeparacion', ...
    'separacionMinima', ...
    'separacionMedia', ...
    'separacionMaxima', ...
    'indiceDB', ...
    'costeValidado', ...
    'diferenciaCoste', ...
    'convergenciaMonotona', ...
    'reduccionCoste', ...
    'parametros');

%% 15. Figura de medias espectrales reales

fig = figure('Visible', 'off');

bar(1:K, mediasEspectrales, 'grouped');

xlabel('Grupo');
ylabel('Valor medio del indice espectral');

title('Caracterizacion espectral de los grupos - DE');

legend({'NDVI', 'NDWI', 'NDBI'}, ...
    'Location', 'best');

xticks(1:K);
grid on;

saveas(fig, ...
    fullfile(figuresDir, 'de_validacion_espectral.png'));

close(fig);

%% 16. Figura de compactacion por grupo

fig = figure('Visible', 'off');

bar(1:K, compactacion);

hold on;

yline(compactacionGlobal, '--', 'Compactacion global');

xlabel('Grupo');
ylabel('Distancia cuadratica media');

title('Compactacion intragrupo - DE');

xticks(1:K);
grid on;

saveas(fig, ...
    fullfile(figuresDir, 'de_validacion_compactacion.png'));

close(fig);

%% 17. Figura de separacion entre centroides

fig = figure('Visible', 'off');

imagesc(matrizSeparacion);

axis image;

colorbar;

xlabel('Centroide');
ylabel('Centroide');

title('Distancias entre centroides - DE');

xticks(1:K);
yticks(1:K);

saveas(fig, ...
    fullfile(figuresDir, 'de_validacion_separacion.png'));

close(fig);

%% 18. Comparar centroides y medias reales

fig = figure('Visible', 'off');

for d = 1:3

    subplot(1, 3, d);

    bar(1:K, ...
        [centrosOriginales(:, d), mediasEspectrales(:, d)], ...
        'grouped');

    xlabel('Grupo');

    nombres = {'NDVI', 'NDWI', 'NDBI'};

    title(nombres{d});

    xticks(1:K);
    grid on;

    if d == 1

        ylabel('Valor del indice');

        legend({'Centroide DE', 'Media real'}, ...
            'Location', 'best');

    end

end

sgtitle('Centroides optimizados frente a medias reales');

saveas(fig, ...
    fullfile(figuresDir, 'de_validacion_centroides.png'));

close(fig);

%% 19. Resumen final

fprintf('\n============================================\n');
fprintf(' RESUMEN DE VALIDACION FINAL DE DE\n');
fprintf('============================================\n');

fprintf('Semilla seleccionada: %d\n', datos.mejorSemilla);
fprintf('Coste sobre muestra: %.6f\n', datos.mejorCoste);
fprintf('Coste sobre imagen: %.6f\n', costeValidado);
fprintf('Compactacion global: %.6f\n', compactacionGlobal);
fprintf('Separacion media: %.6f\n', separacionMedia);
fprintf('Davies-Bouldin: %.6f\n', indiceDB);
fprintf('Grupos vacios: %d\n', gruposVacios);
fprintf('Pixeles clasificados: %d\n', nPixeles);
fprintf('Pixeles excluidos: %d\n', nPixelesExcluidos);
fprintf('Reduccion del coste: %.2f%%\n', reduccionCoste);

fprintf('\nComprobaciones estructurales y numericas: CORRECTAS\n');

fprintf('\nResultados guardados en:\n%s\n', resultsDir);

fprintf('\nFiguras guardadas en:\n%s\n', figuresDir);

%% FUNCION AUXILIAR: CARGAR BANDA

function banda = cargar_banda_de(carpeta, patron)

    archivos = dir(fullfile(carpeta, patron));

    assert(numel(archivos) == 1, ...
        'Se esperaba exactamente un archivo para: %s', ...
        patron);

    rutaArchivo = fullfile(carpeta, archivos(1).name);

    banda = double(imread(rutaArchivo));

    fprintf('Banda cargada: %s\n', archivos(1).name);

end
