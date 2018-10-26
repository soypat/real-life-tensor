

function[data] = C2_5(var);
% Objetivo: estudiar c?mo cambian,al variar la seccion, las deformaciones y tensiones m?ximas en una estructura si se la
% modela con barras o vigas.
Lx = 1600;
E = 210E6;
beamDofPerNode = 3;
barDofPerNode = 2;
q = -10;
h = var*3; %espesor de la barra o viga
b = h; %altura de la barra o viga
A = b*h;
I = (b*h^3)/12;
nodeCoordinates = [0      0
                   Lx    -1000
                   Lx     0
                  2*Lx  -1500
                  2*Lx    0
                  3*Lx  -1700
                  3*Lx    0
                  4*Lx  -1500
                  4*Lx    0
                  5*Lx    -1000
                  5*Lx     0
                  6*Lx     0 ];
elements = [1   2
            1   3
            2   3
            2   4
            3   4
            3   5
            4   5
            4   6
            5   6 
            5   7
            6   7
            6   8 
            6   9
            7   9
            8   9
            8   10
            8   11
            9   11
            10  11
            10  12
            11  12];
        

nNodes = length(nodeCoordinates);
nElements = length(elements);

beamLoadVector = zeros(nNodes,beamDofPerNode);

barLoadVector = zeros(nNodes,barDofPerNode);

loadedElements = [2 6 10 14 18 21];

beamDOF = reshape(1:1:beamDofPerNode*nNodes,beamDofPerNode,[])';
barDOF = reshape(1:1:barDofPerNode*nNodes,barDofPerNode,[])';

restrictedBeamDOF = zeros(nNodes,beamDofPerNode);
restrictedBeamDOF([1 12],:) = [1 1 0;0 1 0];
restrictedBeamDOF = reshape(restrictedBeamDOF',[],1);
freeBeamDOF = ~restrictedBeamDOF;

restrictedBarDOF = zeros(nNodes,barDofPerNode);
restrictedBarDOF([1 12],:) = [1 1; 0 1];
restrictedBarDOF = reshape(restrictedBarDOF',[],1);
freeBarDOF = ~restrictedBarDOF;

barGlobalK = zeros(nNodes*barDofPerNode);
beamGlobalK = zeros(nNodes*beamDofPerNode);

barD = zeros(nNodes*barDofPerNode,1);
beamD = zeros(nNodes*beamDofPerNode,1);

beamLocalK = zeros(2*beamDofPerNode);
beamGlobalK = zeros(nNodes*beamDofPerNode);
barGlobalK = zeros(nNodes*barDofPerNode);

vectorLongitudes = zeros(nElements,1);
beamTvector = zeros(6,6,nElements);
barTvector = zeros(2,4,nElements);

%% Ensamblado de matrices globales y vectores de carga

for(e = 1:nElements)
    
    elementVector = nodeCoordinates(elements(e,2),:) - nodeCoordinates(elements(e,1),:);
    L = norm(elementVector);
    cosDirectores = elementVector/L;
    vectorLongitudes(e) = L;
    
    barT = [cosDirectores 0 0;0 0 cosDirectores];
    barTvector(:,:,e) = barT;
    beamLambda = [cosDirectores 0; -cosDirectores(2) cosDirectores(1) 0;0 0 1];
    beamT = blkdiag(beamLambda,beamLambda);
    beamTvector(:,:,e) = beamT;
    
    barLocalK = ((E*A)/L)*[1 -1;-1 1];
    
    beamLocalK([2 3 5 6],[2 3 5 6]) = (E*I/(L^3))*[  12    6*L   -12  6*L
                                                     6*L   4*L^2 -6*L 2*L^2
                                                    -12   -6*L    12  -6*L
                                                     6*L   2*L^2 -6*L 4*L^2];
    
    beamLocalK([1 4],[1 4]) = barLocalK;
    
    beamKrotada = beamT'*beamLocalK*beamT;
    barKrotada = barT'*barLocalK*barT;
    
    dofBars = [barDOF(elements(e,1),:) barDOF(elements(e,2),:)];
    dofBeams = [beamDOF(elements(e,1),:) beamDOF(elements(e,2),:)];
    
    beamGlobalK(dofBeams,dofBeams) = beamGlobalK(dofBeams,dofBeams) + beamKrotada;
    barGlobalK(dofBars,dofBars) = barGlobalK(dofBars,dofBars) + barKrotada;
    
    element2Load = loadedElements(e == loadedElements);
    
    barLoadVector(elements(element2Load,1),2) = barLoadVector(elements(element2Load,1),2) + q*L/2;
    barLoadVector(elements(element2Load,2),2) = barLoadVector(elements(element2Load,2),2) + q*L/2;
    
    
    beamLoadVector(elements(element2Load,1),2) = beamLoadVector(elements(element2Load,1),2) + q*L/2;
    beamLoadVector(elements(element2Load,1),3) = beamLoadVector(elements(element2Load,1),3) + q*L^2/12;
    
    beamLoadVector(elements(element2Load,2),2) = beamLoadVector(elements(element2Load,2),2) + q*L/2;
    beamLoadVector(elements(element2Load,2),3) = beamLoadVector(elements(element2Load,2),3) - q*L^2/12;
end

beamLoadVector = reshape(beamLoadVector',[],1);
barLoadVector = reshape(barLoadVector',[],1);
%% Reducci?n de Matrices
reducedBeamK = beamGlobalK(freeBeamDOF,freeBeamDOF);
reducedBeamLoadVector = beamLoadVector(freeBeamDOF);


reducedBarK = barGlobalK(freeBarDOF,freeBarDOF);
reducedBarLoadVector = barLoadVector(freeBarDOF);

%% C?lculo de desplazamientos

reducedBeamD = reducedBeamK\reducedBeamLoadVector;

reducedBarD = reducedBarK\reducedBarLoadVector;

beamD(freeBeamDOF) = reducedBeamD;
barD(freeBarDOF) = reducedBarD;
%drawBars(nodeCoordinates,elements,barD)
%drawBars(nodeCoordinates,elements,beamD(
barD = reshape(barD,2,[])';
beamD = reshape(beamD,3,[])';







%% Tensiones

beamDivisions = 100;
sigmaBar = zeros(nElements,1);
sigmaNormalBeam = zeros(nElements,1);
sigmaFlexBeam = zeros(nElements,beamDivisions);
sigmaMaxBeam = zeros(nElements,1);



for(e = 1:nElements)
    
    L = vectorLongitudes(e);
    B = [-(1/L) 1/L];
    
    % Barras
   
    barDisplacements = [barD(elements(e,1),:) barD(elements(e,2),:)]';
    localBarDisplacements = barTvector(:,:,e)*barDisplacements;
    sigmaBar(e) = E*B*localBarDisplacements;
    
    %Vigas
    
    X = linspace(0,vectorLongitudes(e),beamDivisions);
    
    beamDisplacements = [beamD(elements(e,1),:) beamD(elements(e,2),:)]';
    
    localBeamDisplacements = beamTvector(:,:,e)*beamDisplacements;
    localBeamAx = localBeamDisplacements([1  4]);
    sigmaNormalBeam(e) = E*B*localBeamAx;
    v1 = localBeamDisplacements(2);
    theta1 = localBeamDisplacements(3);
    v2 = localBeamDisplacements(5);
    theta2 = localBeamDisplacements(6);
    tensionFlex = @(x) (E*h/2)*( ((-6/L^2)+((12*x)/L^3))*v1 + ((-4/L)+((6*x)/L^2))*theta1 + (6/L^2 - ((12*x)/L^3))*v2 + (-2/L + ((6*x)/L^2))*theta2 );
    sigmaFlexBeam(e,:) = tensionFlex(X);
    sigmaMaxBeam(e) = sigmaNormalBeam(e) + sign(sigmaNormalBeam(e))*max(abs(sigmaFlexBeam(e,:)));  
    
end

data.sigma.beam = sigmaMaxBeam;
data.sigma.bar = sigmaBar;
data.displacements.original = nodeCoordinates;
data.displacements.bars = barD;
data.displacements.beams = beamD;
data.geometry.length = vectorLongitudes;
data.geometry.height = h;


end
 