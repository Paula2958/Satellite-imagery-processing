
%% DE_MEJOR_RESULTADO
% Recupera la mejor solucion del analisis de robustez de DE.
% Aplica los mejores centroides a toda la imagen Sentinel-2.
% No vuelve a ejecutar el algoritmo de optimizacion.

clearvars;
clc;
close all;

%% 1. Definir rutas

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));

rawDataDir = fullfile(projectRoot, 'data', 'raw');
resultsDir = fullfile(projectRoot, 'results');
figuresDir = fullfile(resultsDir, 'figures', 'DE');

addpath(scriptDir);

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

fprintf('============================================\n');
fprintf(' MEJOR RESULTADO DE EVOLUCION DIFERENCIAL\n');
fprintf('============================================\n');

%% 2. Cargar resultados de robustez

archivoRobustez = fullfile(resultsDir, 'robustez_de.mat');

assert(exist(archivoRobustez, 'file') == 2, ...
    'Ejecuta primero robustez_de.m.');

datos = load(archivoRobustez);

parametros = datos.parametros;

%% 3. Seleccionar la mejor ejecucion

[mejorCoste, mejorEjecucion] = min(datos.costes);

mejorSemilla = datos.semillas(mejorEjecucion);
mejoresCentros = datos.centroides{mejorEjecucion};
mejorHistorial = datos.historiales{mejorEjecucion};
mejorTiempo = datos.tiempos(mejorEjecucion);

fprintf('\nMejor ejecucion: %d\n', mejorEjecucion);
fprintf('Semilla: %d\n', mejorSemilla);
fprintf('Mejor coste: %.6f\n', mejorCoste);
fprintf('Tiempo de optimizacion: %.4f segundos\n', ...
    mejorTiempo);

%% 4. Cargar bandas Sentinel-2

B02 = cargar_banda_de(rawDataDir, '*_B02_(Raw).tiff');
B03 = cargar_banda_de(rawDataDir, '*_B03_(Raw).tiff');
B04 = cargar_banda_de(rawDataDir, '*_B04_(Raw).tiff');
B08 = cargar_banda_de(rawDataDir, '*_B08_(Raw).tiff');
B11 = cargar_banda_de(rawDataDir, '*_B11_(Raw).tiff');

assert(isequal(size(B02), size(B03), size(B04), ...
    size(B08), size(B11)), ...
    'Las bandas no tienen las mismas dimensiones.');

%% 5. Crear mascara de pixeles validos

mascara = isfinite(B02) & isfinite(B03) & ...
          isfinite(B04) & isfinite(B08) & ...
          isfinite(B11) & ...
          B02 > 0 & B03 > 0 & B04 > 0 & ...
          B08 > 0 & B11 > 0;

%% 6. Calcular indices espectrales

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

%% 7. Construir matriz de datos

pixelesMascara = find(mascara);

Xtotal = [NDVI(mascara), ...
          NDWI(mascara), ...
          NDBI(mascara)];

filasFinitas = all(isfinite(Xtotal), 2);

X = Xtotal(filasFinitas, :);
pixelesValidos = pixelesMascara(filasFinitas);

assert(~isempty(X), ...
    'No hay pixeles validos.');

fprintf('\nPixeles validos: %d\n', size(X, 1));

%% 8. Normalizar datos

media = mean(X, 1);
desv = std(X, 0, 1);

desv(desv == 0) = 1;

Xnormalizado = (X - media) ./ desv;

%% 9. Verificar la muestra original

Xm = Xnormalizado(datos.idx, :);

assert(isequal(size(Xm), size(datos.Xm)), ...
    'Las dimensiones de la muestra no coinciden.');

diferenciaMuestra = max(abs(Xm(:) - datos.Xm(:)));

assert(diferenciaMuestra < 1e-10, ...
    'La muestra no coincide con la utilizada en robustez.');

fprintf('Muestra original verificada correctamente.\n');

%% 10. Validar el coste de la mejor solucion

D_muestra = distancias_de(Xm, mejoresCentros);

[dminMuestra, etiquetasMuestra] = min(D_muestra, [], 2);

conteoMuestra = accumarray(etiquetasMuestra, 1, ...
    [parametros.K, 1]);

gruposVaciosMuestra = sum(conteoMuestra == 0);

costeVerificado = mean(dminMuestra) + ...
    parametros.penalizacion * gruposVaciosMuestra;

assert(abs(costeVerificado - mejorCoste) < 1e-10, ...
    'El coste calculado no coincide con el coste guardado.');

fprintf('Coste verificado: %.6f\n', costeVerificado);
fprintf('Grupos vacios en la muestra: %d\n', ...
    gruposVaciosMuestra);

%% 11. Clasificar todos los pixeles validos

fprintf('\nClasificando la imagen completa...\n');

nPixeles = size(Xnormalizado, 1);

etiquetas = zeros(nPixeles, 1);
distanciasMinimas = zeros(nPixeles, 1);

tamanoBloque = 10000;

for inicio = 1:tamanoBloque:nPixeles

    fin = min(inicio + tamanoBloque - 1, nPixeles);

    Xbloque = Xnormalizado(inicio:fin, :);

    D = distancias_de(Xbloque, mejoresCentros);

    [distanciasMinimas(inicio:fin), ...
        etiquetas(inicio:fin)] = min(D, [], 2);

end

%% 12. Crear mapa de segmentacion

mapaClases = zeros(size(B02), 'uint8');

mapaClases(pixelesValidos) = uint8(etiquetas);

conteoTotal = accumarray(etiquetas, 1, ...
    [parametros.K, 1]);

porcentajes = 100 * conteoTotal / nPixeles;

gruposVaciosTotal = sum(conteoTotal == 0);

costeTotal = mean(distanciasMinimas) + ...
    parametros.penalizacion * gruposVaciosTotal;

%% 13. Recuperar centroides en escala original

centrosOriginales = mejoresCentros .* desv + media;

%% 14. Crear tabla resumen de grupos

tablaResumen = table( ...
    (1:parametros.K)', ...
    conteoTotal, ...
    porcentajes, ...
    centrosOriginales(:, 1), ...
    centrosOriginales(:, 2), ...
    centrosOriginales(:, 3), ...
    'VariableNames', { ...
        'Clase', ...
        'Pixeles', ...
        'Porcentaje', ...
        'NDVI_Centroide', ...
        'NDWI_Centroide', ...
        'NDBI_Centroide'});

fprintf('\nDistribucion de grupos:\n');
disp(tablaResumen);

%% 15. Guardar clasificacion por pixel

[fila, columna] = ind2sub(size(B02), pixelesValidos);

tablaClases = table( ...
    fila, columna, etiquetas, ...
    'VariableNames', {'Fila', 'Columna', 'Clase'});

writetable(tablaClases, ...
    fullfile(resultsDir, 'de_mejor_clases.csv'));

writetable(tablaResumen, ...
    fullfile(resultsDir, 'de_mejor_resumen.csv'));

%% 16. Guardar resultados completos

save(fullfile(resultsDir, 'de_mejor_resultado.mat'), ...
    'mejorEjecucion', ...
    'mejorSemilla', ...
    'mejorCoste', ...
    'mejorTiempo', ...
    'mejoresCentros', ...
    'centrosOriginales', ...
    'mejorHistorial', ...
    'costeVerificado', ...
    'costeTotal', ...
    'gruposVaciosMuestra', ...
    'gruposVaciosTotal', ...
    'conteoMuestra', ...
    'conteoTotal', ...
    'porcentajes', ...
    'mapaClases', ...
    'etiquetas', ...
    'pixelesValidos', ...
    'parametros', ...
    'media', ...
    'desv', ...
    'tablaResumen');

%% 17. Generar mapa de la mejor segmentacion

fig = figure('Visible', 'off');

imagesc(mapaClases);

axis image;
axis ij;

% Negro para pixeles no validos.
% Un color diferente para cada grupo.

colormap([0 0 0; lines(parametros.K)]);

caxis([-0.5 parametros.K + 0.5]);

cb = colorbar;
cb.Ticks = 0:parametros.K;
cb.TickLabels = [{'No valido'}, ...
    arrayfun(@(k) sprintf('Grupo %d', k), ...
    1:parametros.K, 'UniformOutput', false)];

xlabel('Columna');
ylabel('Fila');

title(sprintf('Mejor segmentacion DE - Semilla %d', ...
    mejorSemilla));

saveas(fig, fullfile(figuresDir, 'de_mejor_mapa.png'));

close(fig);

%% 18. Generar distribucion de grupos

fig = figure('Visible', 'off');

b = bar(1:parametros.K, conteoTotal);

b.FaceColor = 'flat';
b.CData = lines(parametros.K);

xlabel('Grupo');
ylabel('Numero de pixeles');

title('Distribucion de grupos - Mejor solucion DE');

xticks(1:parametros.K);
grid on;

saveas(fig, ...
    fullfile(figuresDir, 'de_mejor_distribucion.png'));

close(fig);

%% 19. Resumen final

fprintf('\n============================================\n');
fprintf(' RESUMEN DE LA MEJOR SOLUCION DE DE\n');
fprintf('============================================\n');

fprintf('Mejor ejecucion: %d\n', mejorEjecucion);
fprintf('Mejor semilla: %d\n', mejorSemilla);
fprintf('Coste de optimizacion: %.6f\n', mejorCoste);
fprintf('Coste verificado: %.6f\n', costeVerificado);
fprintf('Coste sobre imagen completa: %.6f\n', costeTotal);
fprintf('Tiempo de optimizacion: %.4f segundos\n', mejorTiempo);
fprintf('Pixeles clasificados: %d\n', nPixeles);
fprintf('Grupos vacios en muestra: %d\n', gruposVaciosMuestra);
fprintf('Grupos vacios en imagen: %d\n', gruposVaciosTotal);

fprintf('\nResultados guardados en:\n%s\n', resultsDir);
fprintf('\nFiguras guardadas en:\n%s\n', figuresDir);

%% FUNCION AUXILIAR: CARGAR BANDA

function banda = cargar_banda_de(carpeta, patron)

    archivos = dir(fullfile(carpeta, patron));

    assert(numel(archivos) == 1, ...
        'Se esperaba exactamente un archivo para: %s', ...
        patron);

    rutaArchivo = fullfile(carpeta, archivos(1).name);

    banda = double(imread(rutaArchivo));

    fprintf('Banda cargada: %s\n', archivos(1).name);

end
