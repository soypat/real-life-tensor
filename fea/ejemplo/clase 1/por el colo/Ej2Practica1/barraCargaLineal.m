function[state]=barraCargaLineal(nElementos,loadFunction)



%%Datos geom?tricos y del material
%nElementos=100; %borrar
area = 2; % en in2
totalLength = 60;  % en in
E = 30E6;  % psi

%% Coordenadas, elementos y funcion de carga

x1 = 1:1:nElementos;
x2 = 2:1:nElementos+1;

dofPerNode = 1;

elementos = [x1' x2'];

nNodos = nElementos+1;

nodos = linspace(0,totalLength,nNodos);

Px =loadFunction; % Carga variable dependiente de x
%Px = @(x) -10*x;
sigmaVector=sparse(nElementos);


%% Construccion de la matriz global y el vector de cargas


kGenerica =[1 -1;-1 1];
kGlobal = sparse(dofPerNode*nNodos,dofPerNode*nNodos);
loadVector = sparse(nNodos,1);


for(i = 1:nElementos)
    
    longitudElemento=nodos(i+1)-nodos(i);
    k = ((area*E)/longitudElemento);
    kLocal =k*kGenerica;
    
    kGlobal([i,i+1],[i,i+1]) = kGlobal([i,i+1],[i,i+1]) + kLocal;
    
    
    q = (Px(nodos(i))+ Px(nodos(i+1)))/2;  % Carga promedio sobre el elemento
    loadVector(i) = loadVector(i) + q/2;    % Divido la carga sobre el elemento entre sus nodos
    loadVector(i+1) = loadVector(i+1) + q/2;  
end


%% Desplazamientos

kReducida = kGlobal(1:end-1,1:end-1);
loadVectorReducido = loadVector(1:end-1);
D = sparse(nNodos,1);
DReducida = sparse(nNodos-1,1);
DReducida = kReducida\loadVectorReducido;
u = DReducida(1);

%% Diagrama de las cargas sobre la barra
hold on

plot(nodos,loadVector)

xlabel('Posicion')
ylabel('Carga')


axis([0 totalLength -800 0])
grid
hold off

%% Tensiones

D(2:end) = DReducida;
loadVector = kGlobal*D;
generalB = [-1 1];
sigma=0; 

for(i = 1:nElementos)
    
    longitudElemento=nodos(i+1)-nodos(i);
    sigmaVector(i) = (1/longitudElemento)*generalB*[D(i);D(i+1)];

end

sigma = sum(sigmaVector(2:end));


state=[u sigma];
end

    
    
    
 





