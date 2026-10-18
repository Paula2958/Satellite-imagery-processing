
%% OH - Sentinel-2 Telde
% Clasificacion definitiva usando la mejor semilla SA

clear; clc; close all;

%% 1. Cargar datos y mejor solucion
cargar_visualizar;
load('results/robustez_sa.mat', ...
    'mejoresCentrosGlobal','mejorCosteGlobal', ...
    'mejorSemilla','media','desv','K','Xm');

Xtotal = [NDVI(mascara), NDWI(mascara), NDBI(mascara)];
validos = all(isfinite(Xtotal),2);

%% 2. Clasificar todos los pixeles validos
etiquetas = zeros(size(Xtotal,1),1);
posValidas = find(validos);

bloque = 10000;

for inicio = 1:bloque:numel(posValidas)

    fin = min(inicio+bloque-1,numel(posValidas));
    posiciones = posValidas(inicio:fin);

    datos = (Xtotal(posiciones,:)-media)./desv;

    D = zeros(size(datos,1),K);

    for k = 1:K
        dif = datos-mejoresCentrosGlobal(k,:);
        D(:,k) = sum(dif.^2,2);
    end

    [~,etiquetas(posiciones)] = min(D,[],2);
end

mapa = zeros(size(mascara));
mapa(mascara) = etiquetas;

%% 3. Estadisticas
conteos = zeros(K,1);

for k = 1:K
    conteos(k) = nnz(etiquetas==k);
end

porcentajes = 100*conteos/sum(conteos);

fprintf('\n===== MEJOR RESULTADO SA =====\n');
fprintf('Mejor semilla: %d\n',mejorSemilla);
fprintf('Coste: %.6f\n',mejorCosteGlobal);

for k = 1:K
    fprintf('Clase %d: %.2f %%\n',k,porcentajes(k));
end

%% 4. Crear carpeta
if ~exist('results/figures','dir')
    mkdir('results/figures');
end

%% 5. Mapa definitivo
fig1 = figure('Name','Mejor clasificacion SA');

imagesc(mapa);
axis image off;

colores = [0 0 0; lines(K)];
colormap(colores);
clim([-0.5 K+0.5]);

cb = colorbar;
cb.Ticks = 0:K;
cb.TickLabels = [{'Sin datos'}, ...
    arrayfun(@(k)sprintf('Clase %d',k), ...
    1:K,'UniformOutput',false)];

title('Clasificacion espectral - Mejor SA');

exportgraphics(fig1, ...
    'results/figures/sa_mejor_mapa.png', ...
    'Resolution',200);

%% 6. Distribucion de clases
fig2 = figure('Name','Distribucion mejor SA');

bar(1:K,porcentajes);
xlabel('Clase espectral');
ylabel('Pixeles validos (%)');
title('Distribucion de clases - Mejor SA');
xticks(1:K);
ylim([0 100]);
grid on;

exportgraphics(fig2, ...
    'results/figures/sa_mejor_distribucion.png', ...
    'Resolution',200);

%% 7. Guardar datos
save('results/sa_mejor_resultado.mat', ...
    'mapa','mejoresCentrosGlobal', ...
    'mejorCosteGlobal','mejorSemilla', ...
    'conteos','porcentajes','K');

tablaClases = table((1:K)',conteos,porcentajes, ...
    'VariableNames', ...
    {'Clase','Pixeles','Porcentaje'});

writetable(tablaClases, ...
    'results/sa_mejor_clases.csv');

fprintf('\nClasificacion definitiva guardada.\n');
