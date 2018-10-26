%Grafica deformada de la malla

function plotDef(nodos,desp,color,elType)

% Numeración de los nodos de los elementos:
%  4---7---3
%  |       |
%  8       6
%  |       |
%  1---5---2
%  
% nodos:     matriz de coordenadas nodales.
% D:      Desplazamientos de los dof ordenados.
% color:     string que especifica el color para los bordes de los
%            elementos.

%Esto está hecho para elementos de 8 nodos. Defino los bordes
vNod = [1 5 2
    2 6 3
    3 7 4
    4 8 1];

%Recorro los bordes asignando valores para csi y eta.
%csi es la primer columna y eta la segunda. Donde dice 0 indico que es un
%valor que varía.
borde = [0 -1
    1 0
    0 1
    -1 0];


for iBord = 1:4
    
    if borde(iBord,1) == 0
        csi = linspace(-1,1,1000);
    else
        csi = linspace(borde(iBord,1),borde(iBord,1),1000);
    end
    
    if borde(iBord,2) == 0
        eta = linspace(-1,1,1000);
    else
        eta = linspace(borde(iBord,2),borde(iBord,2),1000);
    end
    
    if strcmp(elType,'Q10') == true || strcmp(elType,'QM10') == true
        %Funciones de forma Q10
        N5 = (1-csi.^2) .* (1-eta) / 2;
        N6 = (1+csi) .* (1-eta.^2) / 2;
        N7 = (1-csi.^2) .* (1+eta) / 2;
        N8 = (1-csi) .* (1-eta.^2) / 2;
        N1 = (1 - csi).*(1 - eta) / 4 - (N5+N8)/2;
        N2 = (1 + csi).*(1 - eta) / 4 - (N5+N6)/2;
        N3 = (1 + csi).*(1 + eta) / 4 - (N6+N7)/2;
        N4 = (1 - csi).*(1 + eta) / 4 - (N7+N8)/2;
        N9 = csi.*(1 + csi).*(1 - csi)/(2*sqrt(3)/9); %revisar esto! ver por qué constante multiplica
        N10 = eta.*(1 + eta).*(1 - eta)/(2*sqrt(3)/9);
        N = [N1; N2; N3; N4; N5; N6; N7; N8; N9; N10]';
        xOrig = zeros(10,1);
        yOrig = zeros(10,1);
    elseif strcmp(elType,'Q8') == true
        %Funciones de forma Q8
        N5 = (1-csi.^2) .* (1-eta) / 2;
        N6 = (1+csi) .* (1-eta.^2) / 2;
        N7 = (1-csi.^2) .* (1+eta) / 2;
        N8 = (1-csi) .* (1-eta.^2) / 2;
        N1 = (1 - csi).*(1 - eta) / 4 - (N5+N8)/2;
        N2 = (1 + csi).*(1 - eta) / 4 - (N5+N6)/2;
        N3 = (1 + csi).*(1 + eta) / 4 - (N6+N7)/2;
        N4 = (1 - csi).*(1 + eta) / 4 - (N7+N8)/2;
        N = [N1; N2; N3; N4; N5; N6; N7; N8]';
        xOrig = zeros(8,1);
        yOrig = zeros(8,1);
    end
    
    Dreshape = reshape(desp',2,[])';
    
    xOrig(1:8,1) = nodos(:,1);
    xDisp = Dreshape(:,1);
    
    yOrig(1:8,1) = nodos(:,2);
    yDisp = Dreshape(:,2);
    
    x = N*(xDisp + xOrig);
    y = N*(yDisp + yOrig);
    hold on
    plot(x,y,color)
%     plot(x,y,color,'LineWidth',1.5)
%     axis square
end
    
