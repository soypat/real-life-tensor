%% Comparación entre estrucura formada por barras o vigas
clear ; close all ; clc
%% Dimesiones
n = 100;
Lado = linspace(1,1732,n);
%% Cosientes
VectorBarras = zeros(n,3); % [Desplazamientos TensionesMax ElementoMasTensionado]
VectorVigas = zeros(n,3);  % [Desplazamientos TensionesMax ElementoMasTensionado]
CocienteDesplazamientos = zeros(n,1);
CocienteTensiones = zeros(n,1);
for i = 1:n
    [desp tens elem] = PuenteBarras(Lado(i));
    VectorBarras(i,:) = [desp tens elem];
    [desp tens elem] = PuenteVigas(Lado(i));
    VectorVigas(i,:) = [desp tens elem];
    CocienteDesplazamientos(i) = VectorVigas(i,1)/VectorBarras(i,1);
    CocienteTensiones(i) = VectorBarras(i,2)/VectorVigas(i,2);
end