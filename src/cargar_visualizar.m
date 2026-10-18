%% OH - Sentinel-2 Telde: carga y visualizacion inicial
% Colocar este script junto a los seis TIFF descargados de Copernicus.
clear; clc; close all;
archivos = dir('*_B02_(Raw).tiff');
assert(numel(archivos)==1,'No se encuentra un unico TIFF B02 en la carpeta actual.');
prefijo = extractBefore(string(archivos(1).name),'B02_(Raw).tiff');
bandas = ["B02","B03","B04","B08","B11"];
datos = struct();
for k = 1:numel(bandas)
    nombre = char(prefijo + bandas(k) + "_(Raw).tiff");
    assert(isfile(nombre), 'Falta el archivo %s', nombre);
    datos.(char(bandas(k))) = double(imread(nombre));
end
B02=datos.B02; B03=datos.B03; B04=datos.B04;
B08=datos.B08; B11=datos.B11;
mascara = isfinite(B02)&isfinite(B03)&isfinite(B04)& ...
          isfinite(B08)&isfinite(B11)& ...
          (B02>0)&(B03>0)&(B04>0)&(B08>0)&(B11>0);
fprintf('Dimensiones: %d x %d\n',size(B02,1),size(B02,2));
fprintf('Pixeles validos (mascara conservadora): %.2f %%\n',100*nnz(mascara)/numel(mascara));
% Escalado por percentiles en pixeles validos para evitar imagen lavada.
rgb = cat(3,B04,B03,B02);
rgbVis = zeros(size(rgb));
for c=1:3
    canal = rgb(:,:,c);
    v = canal(mascara);
    assert(~isempty(v),'No hay pixeles validos: revisar las bandas.');
    limites = prctile(v,[2 98]);
    rgbVis(:,:,c)=min(max((canal-limites(1))/max(limites(2)-limites(1),eps),0),1);
end
rgbVis(repmat(~mascara,1,1,3))=0;
figure('Name','RGB Sentinel-2');
image(rgbVis); axis image off; title('RGB (B04-B03-B02) - 28/09/2026');
figure('Name','Mascara de datos');
imagesc(mascara); axis image off; colormap(gray); title('Mascara de pixeles validos');
% Indices espectrales: evitar divisiones por cero y excluir zonas invalidas.
NDVI = (B08-B04)./(B08+B04+eps);
NDWI = (B03-B08)./(B03+B08+eps); % McFeeters
NDBI = (B11-B08)./(B11+B08+eps);
indices = {NDVI,NDWI,NDBI}; nombres = {'NDVI','NDWI','NDBI'};
for k=1:3
    I=indices{k}; I(~mascara)=NaN;
    figure('Name',nombres{k});
    imagesc(I,[-1 1]); axis image off; colorbar;
    title(nombres{k} + " - 28/09/2026");
end
