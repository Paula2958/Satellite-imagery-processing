
%% OPTIMIZAR_DE
% Ejecuta Evolucion Diferencial para segmentar la imagen Sentinel-2.
% Variante: DE/rand/1/bin.
%
% La funcion objetivo y el preprocesamiento siguen el planteamiento
% utilizado en optimizar_pso.m.

clearvars;
clc;
close all;

%% 1. Definir rutas

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));

rawDataDir = fullfile(projectRoot, 'data', 'raw');
resultsDir = fullfile(projectRoot, 'results');
figuresDir = fullfile(resultsDir, 'figures', 'DE');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

addpath(scriptDir);

fprintf('============================================\n');
fprintf(' OPTIMIZACION MEDIANTE EVOLUCION DIFERENCIAL\n');
fprintf('============================================\n');

%% 2. Cargar bandas Sentinel-2

B02 = cargar_banda_de(rawDataDir, '*_B02_(Raw).tiff');
B03 = cargar_banda_de(rawDataDir, '*_B03_(Raw).tiff');
B04 = cargar_banda_de(rawDataDir, '*_B04_(Raw).tiff');
B08 = cargar_banda_de(rawDataDir, '*_B08_(Raw).tiff');
B11 = cargar_banda_de(rawDataDir, '*_B11_(Raw).tiff');

assert(isequal(size(B02), size(B03), size(B04), ...
    size(B08), size(B11)), ...
    'Las bandas Sentinel-2 no tienen las mismas dimensiones.');

fprintf('Dimensiones de imagen: %d x %d\n', ...
    size(B02, 1), size(B02, 2));

%% 3. Crear mascara de pixeles validos

mascara = isfinite(B02) & isfinite(B03) & ...
          isfinite(B04) & isfinite(B08) & ...
          isfinite(B11) & ...
          B02 > 0 & B03 > 0 & B04 > 0 & ...
          B08 > 0 & B11 > 0;

fprintf('Pixeles validos en mascara: %d (%.2f%%)\n', ...
    nnz(mascara), 100 * nnz(mascara) / numel(mascara));

%% 4. Calcular indices espectrales

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

%% 5. Construir matriz de datos

pixelesMascara = find(mascara);

Xtotal = [NDVI(mascara), ...
          NDWI(mascara), ...
          NDBI(mascara)];

filasFinitas = all(isfinite(Xtotal), 2);

X = Xtotal(filasFinitas, :);
pixelesValidos = pixelesMascara(filasFinitas);

fprintf('Pixeles utilizados: %d\n', size(X, 1));

assert(~isempty(X), ...
    'No hay pixeles validos para ejecutar la optimizacion.');

%% 6. Normalizar indices mediante z-score

media = mean(X, 1);
desv = std(X, 0, 1);

desv(desv == 0) = 1;

Xnormalizado = (X - media) ./ desv;

%% 7. Seleccionar muestra reproducible

rng(42);

nMuestra = min(3000, size(Xnormalizado, 1));
idx = randperm(size(Xnormalizado, 1), nMuestra);

Xm = Xnormalizado(idx, :);

fprintf('Tamano de muestra para DE: %d pixeles\n', nMuestra);

%% 8. Configurar Evolucion Diferencial

parametros.K = 4;
parametros.NP = 30;
parametros.F = 0.8;
parametros.CR = 0.9;
parametros.maxIter = 200;
parametros.penalizacion = 2;

fprintf('\nParametros DE:\n');
fprintf('K: %d\n', parametros.K);
fprintf('NP: %d\n', parametros.NP);
fprintf('F: %.2f\n', parametros.F);
fprintf('CR: %.2f\n', parametros.CR);
fprintf('Generaciones: %d\n', parametros.maxIter);
fprintf('Penalizacion: %.2f\n', parametros.penalizacion);

%% 9. Ejecutar Evolucion Diferencial

fprintf('\nIniciando optimizacion DE...\n');

resultado = de_core(Xm, parametros);

mejoresCentros = resultado.mejoresCentros;
mejorCosteGlobal = resultado.mejorCosteGlobal;
historial = resultado.historial;
tiempo = resultado.tiempo;

fprintf('\nOptimizacion finalizada.\n');
fprintf('Mejor coste: %.6f\n', mejorCosteGlobal);
fprintf('Tiempo de ejecucion: %.2f segundos\n', tiempo);

%% 10. Evaluar la solucion sobre la muestra

D_muestra = distancias_de(Xm, mejoresCentros);
[~, etiquetasMuestra] = min(D_muestra, [], 2);

conteoMuestra = accumarray(etiquetasMuestra, 1, ...
    [parametros.K, 1]);

gruposVaciosMuestra = sum(conteoMuestra == 0);

fprintf('\nDistribucion de la muestra:\n');

for k = 1:parametros.K
    fprintf('Grupo %d: %d pixeles\n', k, conteoMuestra(k));
end

fprintf('Grupos vacios en la muestra: %d\n', ...
    gruposVaciosMuestra);

%% 11. Clasificar todos los pixeles validos por bloques

fprintf('\nClasificando la imagen completa...\n');

nPixeles = size(Xnormalizado, 1);
etiquetas = zeros(nPixeles, 1);

tamanoBloque = 10000;

for inicio = 1:tamanoBloque:nPixeles

    fin = min(inicio + tamanoBloque - 1, nPixeles);

    Xbloque = Xnormalizado(inicio:fin, :);

    D = distancias_de(Xbloque, mejoresCentros);

    [~, etiquetas(inicio:fin)] = min(D, [], 2);

end

%% 12. Crear mapa de clases

mapaClases = zeros(size(B02), 'uint8');

mapaClases(pixelesValidos) = uint8(etiquetas);

conteoTotal = accumarray(etiquetas, 1, ...
    [parametros.K, 1]);

gruposVaciosTotal = sum(conteoTotal == 0);

fprintf('\nDistribucion de la imagen completa:\n');

for k = 1:parametros.K
    fprintf('Grupo %d: %d pixeles\n', k, conteoTotal(k));
end

fprintf('Grupos vacios en la imagen: %d\n', ...
    gruposVaciosTotal);

%% 13. Guardar resultados numericos

% Tabla de clasificacion con coordenadas de imagen.
[fila, columna] = ind2sub(size(B02), pixelesValidos);

tablaClases = table( ...
    fila, columna, etiquetas, ...
    'VariableNames', {'Fila', 'Columna', 'Clase'});

writetable(tablaClases, ...
    fullfile(resultsDir, 'de_clases.csv'));

% Historial de convergencia.
generacion = (1:numel(historial))';

tablaConvergencia = table( ...
    generacion, historial, ...
    'VariableNames', {'Generacion', 'MejorCoste'});

writetable(tablaConvergencia, ...
    fullfile(resultsDir, 'de_convergencia.csv'));

% Resumen de resultados.
resumen = struct();
resumen.mejorCoste = mejorCosteGlobal;
resumen.tiempo = tiempo;
resumen.centroides = mejoresCentros;
resumen.conteoClases = conteoTotal;
resumen.gruposVacios = gruposVaciosTotal;
resumen.parametros = parametros;

save(fullfile(resultsDir, 'de_resultados.mat'), ...
    'resultado', 'resumen', 'parametros', ...
    'mejoresCentros', 'mejorCosteGlobal', ...
    'historial', 'tiempo', 'conteoTotal', ...
    'mapaClases', 'media', 'desv');

%% 14. Figura de convergencia

fig = figure('Visible', 'off');

plot(1:numel(historial), historial, 'LineWidth', 1.5);
xlabel('Generacion');
ylabel('Mejor coste');
title('Convergencia de Evolucion Diferencial');
grid on;

saveas(fig, fullfile(figuresDir, 'de_convergencia.png'));
close(fig);

%% 15. Figura de distribucion de clases

fig = figure('Visible', 'off');

bar(1:parametros.K, conteoTotal);
xlabel('Grupo');
ylabel('Numero de pixeles');
title('Distribucion de clases - DE');
xticks(1:parametros.K);
grid on;

saveas(fig, fullfile(figuresDir, 'de_distribucion.png'));
close(fig);

%% 16. Figura del mapa de segmentacion

fig = figure('Visible', 'off');

imagesc(mapaClases);
axis image;
set(gca, 'YDir', 'normal');
colorbar;
caxis([0 parametros.K]);
xlabel('Columna');
ylabel('Fila');
title('Segmentacion Sentinel-2 mediante DE');

saveas(fig, fullfile(figuresDir, 'de_mapa.png'));
close(fig);

%% 17. Resumen final

fprintf('\n============================================\n');
fprintf(' RESUMEN DE EVOLUCION DIFERENCIAL\n');
fprintf('============================================\n');
fprintf('Mejor coste: %.6f\n', mejorCosteGlobal);
fprintf('Tiempo: %.2f segundos\n', tiempo);
fprintf('Pixeles clasificados: %d\n', nPixeles);
fprintf('Grupos vacios: %d\n', gruposVaciosTotal);
fprintf('Resultados guardados en: %s\n', resultsDir);
fprintf('Figuras guardadas en: %s\n', figuresDir);

%% FUNCION AUXILIAR: CARGAR BANDA

function banda = cargar_banda_de(carpeta, patron)

    archivos = dir(fullfile(carpeta, patron));

    assert(numel(archivos) == 1, ...
        'Se esperaba exactamente un archivo para el patron: %s', ...
        patron);

    rutaArchivo = fullfile(carpeta, archivos(1).name);

    banda = double(imread(rutaArchivo));

    fprintf('Banda cargada: %s\n', archivos(1).name);

end
