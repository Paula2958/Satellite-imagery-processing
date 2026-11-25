
%% ROBUSTEZ_DE
% Analisis de robustez de Evolucion Diferencial.
% Se realizan 10 ejecuciones con diferentes semillas.
% Todas las ejecuciones utilizan la misma muestra de pixeles.

clearvars;
clc;
close all;

%% 1. Definir rutas

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));

rawDataDir = fullfile(projectRoot, 'data', 'raw');
resultsDir = fullfile(projectRoot, 'results');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

addpath(scriptDir);

fprintf('============================================\n');
fprintf(' ANALISIS DE ROBUSTEZ - EVOLUCION DIFERENCIAL\n');
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

%% 3. Crear mascara de pixeles validos

mascara = isfinite(B02) & isfinite(B03) & ...
          isfinite(B04) & isfinite(B08) & ...
          isfinite(B11) & ...
          B02 > 0 & B03 > 0 & B04 > 0 & ...
          B08 > 0 & B11 > 0;

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

Xtotal = [NDVI(mascara), ...
          NDWI(mascara), ...
          NDBI(mascara)];

filasFinitas = all(isfinite(Xtotal), 2);
X = Xtotal(filasFinitas, :);

assert(~isempty(X), ...
    'No hay pixeles validos para ejecutar el analisis.');

fprintf('Pixeles validos: %d\n', size(X, 1));

%% 6. Normalizar mediante z-score

media = mean(X, 1);
desv = std(X, 0, 1);

desv(desv == 0) = 1;

Xnormalizado = (X - media) ./ desv;

%% 7. Crear una muestra fija

rng(42);

nMuestra = min(3000, size(Xnormalizado, 1));
idx = randperm(size(Xnormalizado, 1), nMuestra);

Xm = Xnormalizado(idx, :);

fprintf('Tamano de muestra fija: %d pixeles\n', nMuestra);

%% 8. Configurar parametros DE

parametros.K = 4;
parametros.NP = 30;
parametros.F = 0.8;
parametros.CR = 0.9;
parametros.maxIter = 200;
parametros.penalizacion = 2;

%% 9. Configurar ejecuciones

semillas = 1:10;
nEjecuciones = numel(semillas);

costes = zeros(nEjecuciones, 1);
tiempos = zeros(nEjecuciones, 1);
gruposVacios = zeros(nEjecuciones, 1);
iteraciones = parametros.maxIter * ones(nEjecuciones, 1);

centroides = cell(nEjecuciones, 1);
historiales = cell(nEjecuciones, 1);

%% 10. Ejecutar analisis de robustez

fprintf('\nIniciando %d ejecuciones independientes...\n', ...
    nEjecuciones);

for e = 1:nEjecuciones

    fprintf('\n--------------------------------------------\n');
    fprintf('Ejecucion %d/%d | Semilla: %d\n', ...
        e, nEjecuciones, semillas(e));

    rng(semillas(e));

    resultado = de_core(Xm, parametros);

    costes(e) = resultado.mejorCosteGlobal;
    tiempos(e) = resultado.tiempo;
    centroides{e} = resultado.mejoresCentros;
    historiales{e} = resultado.historial;

    % Evaluar grupos vacios en la muestra.
    D = distancias_de(Xm, resultado.mejoresCentros);
    [~, etiquetas] = min(D, [], 2);

    conteo = accumarray(etiquetas, 1, ...
        [parametros.K, 1]);

    gruposVacios(e) = sum(conteo == 0);

    fprintf('Mejor coste: %.6f\n', costes(e));
    fprintf('Tiempo: %.2f segundos\n', tiempos(e));
    fprintf('Grupos vacios: %d\n', gruposVacios(e));

end

%% 11. Calcular estadisticas de robustez

mediaCoste = mean(costes);
desvCoste = std(costes);
minCoste = min(costes);
maxCoste = max(costes);

mediaTiempo = mean(tiempos);
desvTiempo = std(tiempos);

totalGruposVacios = sum(gruposVacios);
ejecucionesConVacios = sum(gruposVacios > 0);

[mejorCoste, mejorEjecucion] = min(costes);

fprintf('\n============================================\n');
fprintf(' ESTADISTICAS DE ROBUSTEZ DE DE\n');
fprintf('============================================\n');

fprintf('Coste medio: %.6f\n', mediaCoste);
fprintf('Desviacion estandar del coste: %.6f\n', desvCoste);
fprintf('Coste minimo: %.6f\n', minCoste);
fprintf('Coste maximo: %.6f\n', maxCoste);

fprintf('Tiempo medio: %.4f segundos\n', mediaTiempo);
fprintf('Desviacion estandar del tiempo: %.4f segundos\n', ...
    desvTiempo);

fprintf('Ejecuciones con grupos vacios: %d/%d\n', ...
    ejecucionesConVacios, nEjecuciones);

fprintf('Total de grupos vacios: %d\n', totalGruposVacios);

fprintf('Mejor ejecucion: %d\n', mejorEjecucion);
fprintf('Semilla de la mejor ejecucion: %d\n', ...
    semillas(mejorEjecucion));
fprintf('Mejor coste obtenido: %.6f\n', mejorCoste);

%% 12. Crear tabla de resultados

tablaRobustez = table( ...
    semillas(:), ...
    costes, ...
    tiempos, ...
    iteraciones, ...
    gruposVacios, ...
    'VariableNames', { ...
        'Semilla', ...
        'Coste', ...
        'Tiempo', ...
        'Iteraciones', ...
        'GruposVacios'});

writetable(tablaRobustez, ...
    fullfile(resultsDir, 'robustez_de.csv'));

%% 13. Guardar resumen y resultados completos

resumen = struct();

resumen.mediaCoste = mediaCoste;
resumen.desviacionCoste = desvCoste;
resumen.minCoste = minCoste;
resumen.maxCoste = maxCoste;

resumen.mediaTiempo = mediaTiempo;
resumen.desviacionTiempo = desvTiempo;

resumen.ejecucionesConVacios = ejecucionesConVacios;
resumen.totalGruposVacios = totalGruposVacios;

resumen.mejorEjecucion = mejorEjecucion;
resumen.mejorCoste = mejorCoste;
resumen.semillaMejorEjecucion = semillas(mejorEjecucion);

save(fullfile(resultsDir, 'robustez_de.mat'), ...
    'tablaRobustez', 'resumen', ...
    'costes', 'tiempos', 'gruposVacios', ...
    'semillas', 'parametros', ...
    'centroides', 'historiales', ...
    'mejorEjecucion', 'Xm', 'idx');

%% 14. Mostrar tabla final

fprintf('\nResultados por ejecucion:\n');
disp(tablaRobustez);

fprintf('\nArchivos guardados en:\n%s\n', resultsDir);

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
