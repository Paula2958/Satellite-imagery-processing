
%% GRAFICAS_ROBUSTEZ_DE
% Genera las figuras del analisis de robustez de DE.

clearvars;
clc;
close all;

%% 1. Definir rutas

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));

resultsDir = fullfile(projectRoot, 'results');
figuresDir = fullfile(resultsDir, 'figures', 'DE');

archivoDatos = fullfile(resultsDir, 'robustez_de.mat');

assert(exist(archivoDatos, 'file') == 2, ...
    'No se encuentra robustez_de.mat. Ejecuta primero robustez_de.m.');

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% 2. Cargar resultados

datos = load(archivoDatos);

costes = datos.costes;
tiempos = datos.tiempos;
semillas = datos.semillas;

mediaCoste = datos.resumen.mediaCoste;
desvCoste = datos.resumen.desviacionCoste;

mediaTiempo = datos.resumen.mediaTiempo;
desvTiempo = datos.resumen.desviacionTiempo;


%% 3. Grafica de distribucion de costes

fig = figure('Visible', 'off');

scatter(semillas, costes, 65, 'filled');
hold on;

yline(mediaCoste, '--', 'Media');

xlabel('Semilla');
ylabel('Coste');
title('Robustez de Evolucion Diferencial: distribucion de costes');

xticks(semillas);
grid on;

legend('Coste por ejecucion', 'Coste medio', ...
    'Location', 'best');

saveas(fig, fullfile(figuresDir, 'robustez_de.png'));
close(fig);


%% 4. Grafica de tiempos por ejecucion

fig = figure('Visible', 'off');

bar(semillas, tiempos);
hold on;

yline(mediaTiempo, '--', 'Tiempo medio', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Semilla');
ylabel('Tiempo (segundos)');
title('Robustez de Evolucion Diferencial: tiempos de ejecucion');
xticks(semillas);
grid on;

saveas(fig, fullfile(figuresDir, 'robustez_de_tiempos.png'));
close(fig);

%% 5. Mostrar resumen

fprintf('============================================\n');
fprintf(' FIGURAS DE ROBUSTEZ DE DE GENERADAS\n');
fprintf('============================================\n');

fprintf('Coste medio: %.6f\n', mediaCoste);
fprintf('Desviacion del coste: %.6f\n', desvCoste);
fprintf('Tiempo medio: %.4f segundos\n', mediaTiempo);
fprintf('Desviacion del tiempo: %.4f segundos\n', desvTiempo);

fprintf('\nFiguras guardadas en:\n%s\n', figuresDir);
