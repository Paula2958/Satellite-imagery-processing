
%% OH - Sentinel-2 Telde
% Validacion final de Simulated Annealing

clear; clc; close all;

%% 1. Cargar resultados
datos = load('results/sa_mejor_resultado.mat');
rob = load('results/robustez_sa.mat');

K = datos.K;
centros = datos.mejoresCentrosGlobal;

%% 2. Recuperar centroides originales
% Los centroides optimizados estaban normalizados
centrosOriginales = centros .* rob.desv + rob.media;

fprintf('\n===== CENTROIDES ESPECTRALES =====\n');

tablaCentros = table( ...
    (1:K)', ...
    centrosOriginales(:,1), ...
    centrosOriginales(:,2), ...
    centrosOriginales(:,3), ...
    datos.porcentajes(:), ...
    'VariableNames', ...
    {'Clase','NDVI','NDWI','NDBI','Porcentaje'});

disp(tablaCentros);

%% 3. Comprobar integridad del mapa
mapa = datos.mapa;

fprintf('\n===== VALIDACION DEL MAPA =====\n');
fprintf('Mejor semilla: %d\n',datos.mejorSemilla);
fprintf('Mejor coste: %.6f\n',datos.mejorCosteGlobal);
fprintf('Numero de clases: %d\n',K);
fprintf('Pixeles clasificados: %d\n',nnz(mapa>0));
fprintf('Clases presentes: %d\n',numel(unique(mapa(mapa>0))));
fprintf('Suma porcentajes: %.4f %%\n',sum(datos.porcentajes));

assert(all(ismember(unique(mapa(:)),0:K)), ...
    'El mapa contiene etiquetas no validas');

assert(all(datos.conteos>0), ...
    'Hay alguna clase vacia');

assert(abs(sum(datos.porcentajes)-100)<0.01, ...
    'Los porcentajes no suman 100');

%% 4. Estadisticas de robustez
costes = rob.costes(:);
tiempos = rob.tiempos(:);

fprintf('\n===== ROBUSTEZ =====\n');
fprintf('Ejecuciones: %d\n',numel(costes));
fprintf('Coste medio: %.6f\n',mean(costes));
fprintf('Desviacion tipica: %.6f\n',std(costes));
fprintf('Mejor coste: %.6f\n',min(costes));
fprintf('Peor coste: %.6f\n',max(costes));
fprintf('Tiempo medio: %.4f s\n',mean(tiempos));

assert(abs(min(costes)-datos.mejorCosteGlobal)<1e-6, ...
    'El coste seleccionado no coincide con robustez');

%% 5. Exportar tabla para la memoria
if ~exist('results','dir')
    mkdir('results');
end

writetable(tablaCentros, ...
    'results/sa_centroides_interpretacion.csv');

fprintf('\nVALIDACION COMPLETADA CORRECTAMENTE\n');
