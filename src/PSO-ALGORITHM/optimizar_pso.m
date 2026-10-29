%% OH - Sentinel-2 Telde
% Script principal de Particle Swarm Optimization (PSO).
% Optimiza cuatro centroides sobre NDVI, NDWI y NDBI.

clear; clc; close all;
rng(42);

%% 1. Rutas del proyecto
algorithmDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(algorithmDir);
projectRoot = fileparts(srcDir);

rawDataDir = fullfile(projectRoot, 'data', 'raw');
figuresDir = fullfile(projectRoot, 'results', 'figures', 'PSO');
metricsDir = fullfile(projectRoot, 'results');

if ~exist(rawDataDir, 'dir')
    error('No se encuentra la carpeta data/raw.');
end
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end
if ~exist(metricsDir, 'dir')
    mkdir(metricsDir);
end

%% 2. Cargar las cinco bandas TIFF
archivos = dir(fullfile(rawDataDir, '*_B02_(Raw).tiff'));

assert(numel(archivos) == 1, ...
    'Se esperaba un unico TIFF B02 en data/raw.');

prefijo = extractBefore(string(archivos(1).name), 'B02_(Raw).tiff');
bandas = ["B02", "B03", "B04", "B08", "B11"];
datos = struct();

for k = 1:numel(bandas)
    nombre = fullfile(rawDataDir, ...
        char(prefijo + bandas(k) + "_(Raw).tiff"));

    assert(isfile(nombre), 'Falta el archivo: %s', nombre);
    datos.(char(bandas(k))) = double(imread(nombre));
end

B02 = datos.B02;
B03 = datos.B03;
B04 = datos.B04;
B08 = datos.B08;
B11 = datos.B11;

assert(isequal(size(B02), size(B03), size(B04), size(B08), size(B11)), ...
    'Las cinco bandas deben tener las mismas dimensiones.');

%% 3. Crear mascara e indices espectrales
mascara = isfinite(B02) & isfinite(B03) & ...
          isfinite(B04) & isfinite(B08) & isfinite(B11) & ...
          B02 > 0 & B03 > 0 & B04 > 0 & ...
          B08 > 0 & B11 > 0;

denNDVI = B08 + B04;
denNDWI = B03 + B08;
denNDBI = B11 + B08;

NDVI = nan(size(B02));
NDWI = nan(size(B02));
NDBI = nan(size(B02));

okNDVI = mascara & isfinite(denNDVI) & denNDVI ~= 0;
okNDWI = mascara & isfinite(denNDWI) & denNDWI ~= 0;
okNDBI = mascara & isfinite(denNDBI) & denNDBI ~= 0;

NDVI(okNDVI) = (B08(okNDVI) - B04(okNDVI)) ./ denNDVI(okNDVI);
NDWI(okNDWI) = (B03(okNDWI) - B08(okNDWI)) ./ denNDWI(okNDWI);
NDBI(okNDBI) = (B11(okNDBI) - B08(okNDBI)) ./ denNDBI(okNDBI);

Xtotal = [NDVI(mascara), NDWI(mascara), NDBI(mascara)];
validos = all(isfinite(Xtotal), 2);

XsinNormalizar = Xtotal(validos, :);
posicionesMascara = find(mascara);
posicionesValidas = posicionesMascara(validos);

K = 4;

assert(size(XsinNormalizar, 1) >= K, ...
    'No hay suficientes pixeles validos para PSO.');

fprintf('Pixeles de la mascara: %d\n', nnz(mascara));
fprintf('Pixeles usados por PSO: %d\n', size(XsinNormalizar, 1));
fprintf('Porcentaje usado: %.2f %%\n', ...
    100 * size(XsinNormalizar, 1) / numel(mascara));

%% 4. Normalizar los indices
media = mean(XsinNormalizar, 1);
desv = std(XsinNormalizar, 0, 1);
desv(desv == 0) = 1;

X = (XsinNormalizar - media) ./ desv;

%% 5. Crear una muestra fija para optimizar
nMuestra = min(3000, size(X, 1));
idx = randperm(size(X, 1), nMuestra);
Xm = X(idx, :);

%% 6. Parametros de PSO
parametros.K = K;
parametros.nParticulas = 30;
parametros.maxIter = 200;
parametros.w = 0.72;
parametros.c1 = 1.49;
parametros.c2 = 1.49;
parametros.penalizacion = 2;

%% 7. Ejecutar PSO
% pso_core.m contiene la logica comun del algoritmo.
resultado = pso_core(Xm, parametros);

mejoresCentros = resultado.mejoresCentros;
mejorCosteGlobal = resultado.mejorCosteGlobal;
historial = resultado.historial;
tiempo = resultado.tiempo;

%% 8. Evaluar la muestra final
D = distancias_pso(Xm, mejoresCentros);
[distanciasMinimas, etiquetasMuestra] = min(D, [], 2);

costeBase = mean(distanciasMinimas);
gruposVacios = K - numel(unique(etiquetasMuestra));

%% 9. Clasificar todos los pixeles validos
% Se clasifican por bloques para reducir el uso de memoria.
etiquetas = zeros(size(X, 1), 1);
bloque = 10000;

for inicio = 1:bloque:size(X, 1)
    fin = min(inicio + bloque - 1, size(X, 1));
    posiciones = inicio:fin;

    D = distancias_pso(X(posiciones, :), mejoresCentros);
    [~, etiquetas(posiciones)] = min(D, [], 2);
end

mapa = zeros(size(mascara));
mapa(posicionesValidas) = etiquetas;

%% 10. Estadisticas
conteos = zeros(K, 1);

for k = 1:K
    conteos(k) = nnz(etiquetas == k);
end

porcentajes = 100 * conteos / sum(conteos);

fprintf('\n===== RESULTADOS PRELIMINARES PSO =====\n');
fprintf('Mejor coste penalizado: %.6f\n', mejorCosteGlobal);
fprintf('Coste sin penalizar: %.6f\n', costeBase);
fprintf('Tiempo: %.4f segundos\n', tiempo);
fprintf('Iteraciones: %d\n', parametros.maxIter);
fprintf('Particulas: %d\n', parametros.nParticulas);
fprintf('Grupos vacios en muestra: %d\n', gruposVacios);

for k = 1:K
    fprintf('Clase %d: %.2f %%\n', k, porcentajes(k));
end

%% 11. Figura de convergencia
fig1 = figure('Name', 'Convergencia PSO');

plot(1:numel(historial), historial, 'LineWidth', 1.5);
xlabel('Iteracion');
ylabel('Mejor coste penalizado');
title('Convergencia del algoritmo PSO');
grid on;

exportgraphics(fig1, ...
    fullfile(figuresDir, 'pso_convergencia.png'), ...
    'Resolution', 200);

%% 12. Figura de distribucion
fig2 = figure('Name', 'Distribucion PSO');

bar(1:K, porcentajes);
xlabel('Clase espectral');
ylabel('Pixeles validos (%)');
title('Distribucion de clases - PSO');
xticks(1:K);
ylim([0 100]);
grid on;

exportgraphics(fig2, ...
    fullfile(figuresDir, 'pso_distribucion.png'), ...
    'Resolution', 200);

%% 13. Figura del mapa
fig3 = figure('Name', 'Mapa PSO');

imagesc(mapa);
axis image off;

colores = [0 0 0; lines(K)];
colormap(colores);
clim([-0.5 K + 0.5]);

cb = colorbar;
cb.Ticks = 0:K;
cb.TickLabels = [{'Sin datos'}, ...
    arrayfun(@(k) sprintf('Clase %d', k), ...
    1:K, 'UniformOutput', false)];

title('Clasificacion espectral - PSO');

exportgraphics(fig3, ...
    fullfile(figuresDir, 'pso_mapa.png'), ...
    'Resolution', 200);

%% 14. Guardar resultados numericos
save(fullfile(metricsDir, 'pso_resultados.mat'), ...
    'mapa', 'mejoresCentros', 'mejorCosteGlobal', ...
    'costeBase', 'tiempo', 'historial', 'conteos', ...
    'porcentajes', 'gruposVacios', 'K', 'media', ...
    'desv', 'idx', 'parametros');

tablaClases = table((1:K)', conteos, porcentajes, ...
    'VariableNames', {'Clase', 'Pixeles', 'Porcentaje'});

writetable(tablaClases, ...
    fullfile(metricsDir, 'pso_clases.csv'));

tablaConvergencia = table((1:numel(historial))', historial, ...
    'VariableNames', {'Iteracion', 'MejorCoste'});

writetable(tablaConvergencia, ...
    fullfile(metricsDir, 'pso_convergencia.csv'));

fprintf('\nFiguras guardadas en results/PSO.\n');
fprintf('Metricas guardadas en results/metrics.\n');
