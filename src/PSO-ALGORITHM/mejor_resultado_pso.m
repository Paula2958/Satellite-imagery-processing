%% OH - Sentinel-2 Telde
% Genera el mapa y la distribucion de la mejor solucion PSO guardada.

clear; clc; close all;

%% 1. Rutas del proyecto
algorithmDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(algorithmDir);
projectRoot = fileparts(srcDir);

metricsDir = fullfile(projectRoot, 'results', 'metrics');
figuresDir = fullfile(projectRoot, 'results', 'PSO');

archivoResultados = fullfile(metricsDir, 'pso_resultados.mat');

if ~isfile(archivoResultados)
    error(['No se encuentra pso_resultados.mat. ', ...
        'Ejecuta primero optimizar_pso.m.']);
end

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% 2. Cargar la mejor solucion guardada
datos = load(archivoResultados, ...
    'mapa', 'mejoresCentros', 'mejorCosteGlobal', ...
    'costeBase', 'tiempo', 'conteos', 'porcentajes', 'K');

mapa = datos.mapa;
mejoresCentros = datos.mejoresCentros;
mejorCosteGlobal = datos.mejorCosteGlobal;
costeBase = datos.costeBase;
tiempo = datos.tiempo;
conteos = datos.conteos;
porcentajes = datos.porcentajes;
K = datos.K;

%% 3. Mostrar estadisticas
fprintf('\n===== MEJOR RESULTADO PSO GUARDADO =====\n');
fprintf('Mejor coste penalizado: %.6f\n', mejorCosteGlobal);
fprintf('Coste sin penalizar: %.6f\n', costeBase);
fprintf('Tiempo de ejecucion: %.4f segundos\n', tiempo);

for k = 1:K
    fprintf('Clase %d: %d pixeles (%.2f %%)\n', ...
        k, conteos(k), porcentajes(k));
end

%% 4. Generar el mapa de la solucion
fig1 = figure('Name', 'Mejor mapa PSO');

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

title('Clasificacion espectral - Mejor PSO');

exportgraphics(fig1, ...
    fullfile(figuresDir, 'pso_mejor_mapa.png'), ...
    'Resolution', 200);

%% 5. Generar la distribucion de clases
fig2 = figure('Name', 'Mejor distribucion PSO');

bar(1:K, porcentajes);
xlabel('Clase espectral');
ylabel('Pixeles validos (%)');
title('Distribucion de clases - Mejor PSO');
xticks(1:K);
ylim([0 100]);
grid on;

exportgraphics(fig2, ...
    fullfile(figuresDir, 'pso_mejor_distribucion.png'), ...
    'Resolution', 200);

%% 6. Guardar resumen de la mejor solucion
tablaClases = table((1:K)', conteos, porcentajes, ...
    'VariableNames', {'Clase', 'Pixeles', 'Porcentaje'});

writetable(tablaClases, ...
    fullfile(metricsDir, 'pso_mejor_clases.csv'));

fprintf('\nMejor mapa guardado en results/PSO.\n');
fprintf('Mejor distribucion guardada en results/PSO.\n');
fprintf('Estadisticas guardadas en results/metrics.\n');
