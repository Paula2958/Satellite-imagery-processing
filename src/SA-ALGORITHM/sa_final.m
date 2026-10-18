
%% OH - Sentinel-2 Telde
% Ejecucion final de Simulated Annealing
% Clasificacion espectral con cuatro centroides

clear; clc; close all;
rng(42);

%% 1. Cargar y preparar datos
cargar_visualizar;

Xtotal = [NDVI(mascara), NDWI(mascara), NDBI(mascara)];
validos = all(isfinite(Xtotal),2);
X = Xtotal(validos,:);

K = 4;
assert(size(X,1)>=K,'Datos insuficientes');

media = mean(X,1);
desv = std(X,0,1);
desv(desv==0) = 1;

X = (X-media)./desv;

% Muestra fija para optimizacion
nMuestra = min(3000,size(X,1));
idx = randperm(size(X,1),nMuestra);
Xm = X(idx,:);

%% 2. Parametros definitivos
T0 = 1.0;
T = T0;
Tmin = 0.001;
alpha = 0.98;
iterPorTemp = 30;
paso = 0.15;
penalizacion = 2;

%% 3. Inicializacion separada
centros = zeros(K,size(Xm,2));
centros(1,:) = Xm(randi(size(Xm,1)),:);

for k = 2:K
    D = distanciasCuadradas(Xm,centros(1:k-1,:));
    [~,pos] = max(min(D,[],2));
    centros(k,:) = Xm(pos,:);
end

costeInicial = costePenalizado(Xm,centros,penalizacion);
costeActual = costeInicial;
mejorCoste = costeActual;
mejoresCentros = centros;

% Historial de convergencia
nNiveles = 0;
tempAux = T0;
while tempAux > Tmin
    nNiveles = nNiveles+1;
    tempAux = alpha*tempAux;
end

maxIter = nNiveles*iterPorTemp;
historial = zeros(maxIter,1);
iteracion = 0;
aceptados = 0;
mejoras = 0;

%% 4. Ejecutar Simulated Annealing
tic;

while T > Tmin

    for j = 1:iterPorTemp
        iteracion = iteracion+1;

        nuevosCentros = centros;
        c = randi(K);

        nuevosCentros(c,:) = nuevosCentros(c,:) ...
            + paso*randn(1,size(Xm,2));

        nuevoCoste = costePenalizado( ...
            Xm,nuevosCentros,penalizacion);

        delta = nuevoCoste-costeActual;

        % Criterio de Metropolis
        if delta <= 0 || rand < exp(-delta/T)
            centros = nuevosCentros;
            costeActual = nuevoCoste;
            aceptados = aceptados+1;
        end

        % Mejor solucion global
        if costeActual < mejorCoste
            mejorCoste = costeActual;
            mejoresCentros = centros;
            mejoras = mejoras+1;
        end

        historial(iteracion) = mejorCoste;
    end

    T = alpha*T;
end

tiempo = toc;
historial = historial(1:iteracion);

%% 5. Coste sin penalizacion
D = distanciasCuadradas(Xm,mejoresCentros);
[costesMinimos,etiquetasMuestra] = min(D,[],2);
costeBase = mean(costesMinimos);

gruposVacios = K-numel(unique(etiquetasMuestra));

%% 6. Clasificar imagen completa
etiquetas = zeros(size(Xtotal,1),1);
posValidas = find(validos);

bloque = 10000;

for inicio = 1:bloque:numel(posValidas)

    fin = min(inicio+bloque-1,numel(posValidas));
    posiciones = posValidas(inicio:fin);

    datos = (Xtotal(posiciones,:)-media)./desv;

    D = distanciasCuadradas(datos,mejoresCentros);
    [~,etiquetas(posiciones)] = min(D,[],2);
end

mapa = zeros(size(mascara));
mapa(mascara) = etiquetas;

%% 7. Estadisticas de clases
conteos = zeros(K,1);

for k = 1:K
    conteos(k) = nnz(etiquetas==k);
end

porcentajes = 100*conteos/sum(conteos);

%% 8. Mostrar resultados
fprintf('\n===== RESULTADOS FINALES SA =====\n');
fprintf('Coste inicial: %.6f\n',costeInicial);
fprintf('Mejor coste penalizado: %.6f\n',mejorCoste);
fprintf('Coste sin penalizar: %.6f\n',costeBase);
fprintf('Iteraciones: %d\n',iteracion);
fprintf('Tiempo: %.4f segundos\n',tiempo);
fprintf('Movimientos aceptados: %d\n',aceptados);
fprintf('Mejoras globales: %d\n',mejoras);
fprintf('Grupos vacios en muestra: %d\n',gruposVacios);

for k = 1:K
    fprintf('Clase %d: %.2f %%\n',k,porcentajes(k));
end

%% 9. Crear carpeta de resultados
if ~exist('results/figures','dir')
    mkdir('results/figures');
end

%% 10. Mapa de clasificacion
fig1 = figure('Name','SA final - Clasificacion');

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

title('Clasificacion espectral - SA final');

exportgraphics(fig1, ...
    'results/figures/sa_final_mapa.png', ...
    'Resolution',200);

%% 11. Grafica de convergencia
fig2 = figure('Name','SA final - Convergencia');

plot(historial,'LineWidth',1.5);
xlabel('Iteracion');
ylabel('Mejor coste penalizado');
title('Convergencia de Simulated Annealing');
grid on;

exportgraphics(fig2, ...
    'results/figures/sa_final_convergencia.png', ...
    'Resolution',200);

%% 12. Distribucion de clases
fig3 = figure('Name','SA final - Distribucion');

bar(1:K,porcentajes);
xlabel('Clase espectral');
ylabel('Pixeles validos (%)');
title('Distribucion de clases - SA final');
xticks(1:K);
ylim([0 100]);
grid on;

exportgraphics(fig3, ...
    'results/figures/sa_final_distribucion.png', ...
    'Resolution',200);

%% 13. Guardar datos completos
save('results/sa_final.mat', ...
    'mapa','mejoresCentros','mejorCoste', ...
    'costeBase','costeInicial','historial', ...
    'conteos','porcentajes','tiempo', ...
    'iteracion','aceptados','mejoras', ...
    'gruposVacios','K','T0','Tmin', ...
    'alpha','paso','penalizacion', ...
    'iterPorTemp','media','desv','idx');

%% 14. CSV para comparacion con PSO y DA
algoritmo = "SA";
tablaResumen = table(algoritmo,K,costeBase, ...
    mejorCoste,tiempo,iteracion, ...
    aceptados,gruposVacios);

writetable(tablaResumen, ...
    'results/sa_final_resumen.csv');

tablaClases = table((1:K)',conteos,porcentajes, ...
    'VariableNames', ...
    {'Clase','Pixeles','Porcentaje'});

writetable(tablaClases, ...
    'results/sa_final_clases.csv');

tablaConvergencia = table( ...
    (1:iteracion)',historial, ...
    'VariableNames',{'Iteracion','MejorCoste'});

writetable(tablaConvergencia, ...
    'results/sa_final_convergencia.csv');

fprintf('\nArchivos finales guardados en results/\n');

%% FUNCIONES AUXILIARES

function D = distanciasCuadradas(X,C)
    D = zeros(size(X,1),size(C,1));

    for k = 1:size(C,1)
        diferencia = X-C(k,:);
        D(:,k) = sum(diferencia.^2,2);
    end
end

function coste = costePenalizado(X,C,penalizacion)
    D = distanciasCuadradas(X,C);

    [distMin,etiquetas] = min(D,[],2);

    costeBase = mean(distMin);

    gruposVacios = size(C,1) ...
        -numel(unique(etiquetas));

    coste = costeBase ...
        +penalizacion*gruposVacios;
end
