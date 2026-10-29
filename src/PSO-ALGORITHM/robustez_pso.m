%% OH - Sentinel-2 Telde
% Analisis de robustez de Particle Swarm Optimization (PSO).
% Compara diez ejecuciones y guarda la mejor solucion encontrada.

clear; clc; close all;

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

%% 4. Normalizar los indices
media = mean(XsinNormalizar, 1);
desv = std(XsinNormalizar, 0, 1);
desv(desv == 0) = 1;

X = (XsinNormalizar - media) ./ desv;

%% 5. Crear una muestra fija para todas las ejecuciones
rng(42);

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

%% 7. Ejecutar PSO con diez semillas
semillas = 1:10;
nEjecuciones = numel(semillas);

costes = zeros(nEjecuciones, 1);
tiempos = zeros(nEjecuciones, 1);

mejorCosteGlobal = inf;
mejorSemilla = NaN;
mejoresCentrosGlobal = [];
mejorHistorial = [];
mejorTiempo = NaN;

fprintf('\n===== ROBUSTEZ PSO =====\n');

for i = 1:nEjecuciones

    rng(semillas(i));

    resultado = pso_core(Xm, parametros);

    costes(i) = resultado.mejorCosteGlobal;
    tiempos(i) = resultado.tiempo;

    fprintf('Semilla %d/%d | coste: %.6f | tiempo: %.4f s\n', ...
        semillas(i), nEjecuciones, costes(i), tiempos(i));

    if costes(i) < mejorCosteGlobal
        mejorCosteGlobal = costes(i);
        mejorSemilla = semillas(i);
        mejoresCentrosGlobal = resultado.mejoresCentros;
        mejorHistorial = resultado.historial;
        mejorTiempo = resultado.tiempo;
    end
end

%% 8. Evaluar la mejor solucion en la muestra
D = distancias_pso(Xm, mejoresCentrosGlobal);
[distanciasMinimas, etiquetasMuestra] = min(D, [], 2);

costeBase = mean(distanciasMinimas);
gruposVacios = K - numel(unique(etiquetasMuestra));

%% 9. Clasificar todos los pixeles validos
etiquetas = zeros(size(X, 1), 1);
bloque = 10000;

for inicio = 1:bloque:size(X, 1)
    fin = min(inicio + bloque - 1, size(X, 1));
    posiciones = inicio:fin;

    D = distancias_pso(X(posiciones, :), mejoresCentrosGlobal);
    [~, etiquetas(posiciones)] = min(D, [], 2);
end

mapa = zeros(size(mascara));
mapa(posicionesValidas) = etiquetas;

%% 10. Estadisticas de la mejor solucion
conteos = zeros(K, 1);

for k = 1:K
    conteos(k) = nnz(etiquetas == k);
end

porcentajes = 100 * conteos / sum(conteos);

fprintf('\n===== MEJOR EJECUCION PSO =====\n');
fprintf('Mejor semilla: %d\n', mejorSemilla);
fprintf('Mejor coste penalizado: %.6f\n', mejorCosteGlobal);
fprintf('Coste sin penalizar: %.6f\n', costeBase);
fprintf('Tiempo de la mejor ejecucion: %.4f segundos\n', mejorTiempo);

%% 11. Figura de robustez
fig1 = figure('Name', 'Robustez PSO');

subplot(2, 1, 1);
plot(semillas, costes, '-o', 'LineWidth', 1.5);
hold on;
plot(mejorSemilla, mejorCosteGlobal, 'p', 'MarkerSize', 12);
xlabel('Semilla aleatoria');
ylabel('Mejor coste penalizado');
title('Coste obtenido en cada ejecucion PSO');
xticks(semillas);
grid on;
hold off;

subplot(2, 1, 2);
bar(semillas, tiempos);
xlabel('Semilla aleatoria');
ylabel('Tiempo (segundos)');
title('Tiempo de cada ejecucion PSO');
xticks(semillas);
grid on;

exportgraphics(fig1, ...
    fullfile(figuresDir, 'pso_robustez.png'), ...
    'Resolution', 200);

%% 12. Guardar tabla de robustez
tablaRobustez = table(semillas(:), costes, tiempos, ...
    'VariableNames', {'Semilla', 'MejorCoste', 'TiempoSegundos'});

writetable(tablaRobustez, ...
    fullfile(metricsDir, 'pso_robustez.csv'));

%% 13. Guardar datos de robustez y mejor solucion
save(fullfile(metricsDir, 'robustez_pso.mat'), ...
    'semillas', 'costes', 'tiempos', ...
    'mejorSemilla', 'mejorCosteGlobal', ...
    'mejoresCentrosGlobal', 'mejorHistorial', ...
    'mejorTiempo', 'costeBase', 'gruposVacios', ...
    'mapa', 'conteos', 'porcentajes', ...
    'media', 'desv', 'K', 'Xm', 'idx', 'parametros');

% Guardar tambien la mejor solucion con el formato
% que utiliza optimizar_pso.m y mejor_resultado_pso.m.
historial = mejorHistorial;
tiempo = mejorTiempo;
mejoresCentros = mejoresCentrosGlobal;

save(fullfile(metricsDir, 'pso_resultados.mat'), ...
    'mapa', 'mejoresCentros', 'mejorCosteGlobal', ...
    'costeBase', 'tiempo', 'historial', 'conteos', ...
    'porcentajes', 'gruposVacios', 'K', 'media', ...
    'desv', 'idx', 'parametros');

tablaClases = table((1:K)', conteos, porcentajes, ...
    'VariableNames', {'Clase', 'Pixeles', 'Porcentaje'});

writetable(tablaClases, ...
    fullfile(metricsDir, 'pso_mejor_clases.csv'));

tablaConvergencia = table((1:numel(historial))', historial, ...
    'VariableNames', {'Iteracion', 'MejorCoste'});

writetable(tablaConvergencia, ...
    fullfile(metricsDir, 'pso_convergencia.csv'));

fprintf('\nMejor semilla: %d\n', mejorSemilla);
fprintf('Resultados de robustez guardados en results/metrics.\n');
fprintf('Figura guardada en results/PSO/pso_robustez.png.\n');
