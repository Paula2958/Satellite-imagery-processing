
%% OH - Sentinel-2 Telde
% Robustez definitiva de Simulated Annealing
% Diez semillas y seleccion de la mejor solucion

clear; clc; close all;
rng(42);

%% 1. Cargar datos
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

% Muestra fija para todas las semillas
nMuestra = min(3000,size(X,1));
idx = randperm(size(X,1),nMuestra);
Xm = X(idx,:);

%% 2. Parametros definitivos
T0 = 1;
Tmin = 0.001;
alpha = 0.98;
iterPorTemp = 30;
paso = 0.15;
penalizacion = 2;

semillas = (1:10)';
nEjecuciones = numel(semillas);

costes = zeros(nEjecuciones,1);
tiempos = zeros(nEjecuciones,1);
iteraciones = zeros(nEjecuciones,1);
vacios = zeros(nEjecuciones,1);
centroides = zeros(K,3,nEjecuciones);

mejorCosteGlobal = inf;
mejorSemilla = NaN;
mejoresCentrosGlobal = [];
mejorHistorial = [];

%% 3. Ejecutar SA con diez semillas
for e = 1:nEjecuciones

    rng(semillas(e));

    centros = zeros(K,size(Xm,2));
    centros(1,:) = Xm(randi(size(Xm,1)),:);

    % Inicializacion por maxima separacion
    for k = 2:K
        D = distanciasSA(Xm,centros(1:k-1,:));
        [~,pos] = max(min(D,[],2));
        centros(k,:) = Xm(pos,:);
    end

    costeActual = costeSA(Xm,centros,penalizacion);
    mejorCoste = costeActual;
    mejoresCentros = centros;

    T = T0;
    iter = 0;
    historial = [];

    tic;

    while T > Tmin

        for j = 1:iterPorTemp

            iter = iter+1;

            nuevos = centros;
            c = randi(K);

            nuevos(c,:) = nuevos(c,:) ...
                + paso*randn(1,size(Xm,2));

            nuevoCoste = costeSA( ...
                Xm,nuevos,penalizacion);

            delta = nuevoCoste-costeActual;

            % Criterio de Metropolis
            if delta <= 0 || rand < exp(-delta/T)
                centros = nuevos;
                costeActual = nuevoCoste;
            end

            if costeActual < mejorCoste
                mejorCoste = costeActual;
                mejoresCentros = centros;
            end

            historial(iter,1) = mejorCoste;
        end

        T = alpha*T;
    end

    tiempos(e) = toc;
    costes(e) = mejorCoste;
    iteraciones(e) = iter;
    centroides(:,:,e) = mejoresCentros;

    D = distanciasSA(Xm,mejoresCentros);
    [~,etiquetas] = min(D,[],2);
    vacios(e) = K-numel(unique(etiquetas));

    % Guardar la mejor ejecucion
    if mejorCoste < mejorCosteGlobal
        mejorCosteGlobal = mejorCoste;
        mejorSemilla = semillas(e);
        mejoresCentrosGlobal = mejoresCentros;
        mejorHistorial = historial;
    end

    fprintf('Semilla %d: coste %.6f | %.3f s\n', ...
        semillas(e),costes(e),tiempos(e));
end

%% 4. Estadisticas generales
costeMedio = mean(costes);
desviacion = std(costes);
mejorCoste = min(costes);
peorCoste = max(costes);
tiempoMedio = mean(tiempos);
ejecucionesVacias = nnz(vacios>0);

fprintf('\n===== ROBUSTEZ DEFINITIVA SA =====\n');
fprintf('Coste medio: %.6f\n',costeMedio);
fprintf('Desviacion tipica: %.6f\n',desviacion);
fprintf('Mejor coste: %.6f\n',mejorCoste);
fprintf('Peor coste: %.6f\n',peorCoste);
fprintf('Tiempo medio: %.4f segundos\n',tiempoMedio);
fprintf('Mejor semilla: %d\n',mejorSemilla);
fprintf('Ejecuciones con grupos vacios: %d\n', ...
    ejecucionesVacias);

%% 5. Crear carpetas
if ~exist('results/figures','dir')
    mkdir('results/figures');
end

%% 6. Guardar CSV de robustez
tabla = table(semillas,costes,tiempos, ...
    iteraciones,vacios, ...
    'VariableNames',{'Semilla','Coste', ...
    'Tiempo','Iteraciones','GruposVacios'});

writetable(tabla,'results/robustez_sa.csv');

%% 7. Guardar datos completos
save('results/robustez_sa.mat', ...
    'tabla','costes','tiempos','semillas', ...
    'iteraciones','vacios','centroides', ...
    'mejoresCentrosGlobal','mejorCosteGlobal', ...
    'mejorSemilla','mejorHistorial', ...
    'Xm','media','desv','idx', ...
    'K','T0','Tmin','alpha', ...
    'iterPorTemp','paso','penalizacion');

%% 8. Grafica de robustez
fig1 = figure('Name','Robustez SA');

plot(semillas,costes,'-o','LineWidth',1.5);
hold on;
yline(costeMedio,'--','Media');

xlabel('Semilla aleatoria');
ylabel('Mejor coste penalizado');
title('Robustez de Simulated Annealing');
xticks(semillas);
grid on;

exportgraphics(fig1, ...
    'results/figures/robustez_sa.png', ...
    'Resolution',200);

%% 9. Grafica de tiempos
fig2 = figure('Name','Tiempos SA');

bar(semillas,tiempos);
xlabel('Semilla aleatoria');
ylabel('Tiempo (s)');
title('Tiempo por ejecucion de SA');
xticks(semillas);
grid on;

exportgraphics(fig2, ...
    'results/figures/robustez_sa_tiempos.png', ...
    'Resolution',200);

%% 10. Convergencia de la mejor semilla
fig3 = figure('Name','Mejor convergencia SA');

plot(mejorHistorial,'LineWidth',1.5);
xlabel('Iteracion');
ylabel('Mejor coste penalizado');
title('Convergencia de la mejor ejecucion SA');
grid on;

exportgraphics(fig3, ...
    'results/figures/robustez_sa_mejor_convergencia.png', ...
    'Resolution',200);

%% 11. Guardar resumen para comparacion
Algoritmo = "SA";
Semillas = nEjecuciones;
CosteMedio = costeMedio;
Desviacion = desviacion;
MejorCoste = mejorCoste;
PeorCoste = peorCoste;
TiempoMedio = tiempoMedio;
MejorSemilla = mejorSemilla;
GruposVacios = ejecucionesVacias;

resumen = table(Algoritmo,Semillas, ...
    CosteMedio,Desviacion,MejorCoste, ...
    PeorCoste,TiempoMedio, ...
    MejorSemilla,GruposVacios);

writetable(resumen, ...
    'results/robustez_sa_resumen.csv');

fprintf('\nArchivos guardados correctamente.\n');

%% FUNCIONES AUXILIARES

function D = distanciasSA(X,C)
    D = zeros(size(X,1),size(C,1));

    for k = 1:size(C,1)
        dif = X-C(k,:);
        D(:,k) = sum(dif.^2,2);
    end
end

function coste = costeSA(X,C,penalizacion)
    D = distanciasSA(X,C);
    [dmin,etiquetas] = min(D,[],2);

    vacios = size(C,1)-numel(unique(etiquetas));
    coste = mean(dmin)+penalizacion*vacios;
end
