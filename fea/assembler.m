%Assembler
cont=true;
nodos=[];

while cont
    stringi=input('ingrese nodo: ','s');
    noditos=sscanf(stringi,'%f');
    if isempty(noditos)
        cont=false;
        break
    end
    nodos=[nodos;noditos(1) noditos(2)]; 
    fprintf('\n')
end
[N placeholder]=size(nodos);
fprintf('N=%0.0f',N);
cont=true;
i=0;
fprintf('\n Datos ingresados: E  A  I\n\nPresione enter para tomar datos anteriores. \n Ingrese nodos en blanco para terminar.\n\n');
Evec=[];
Avec=[];
Ivec=[];
while cont 
    i=i+1;
    inputmessage=sprintf('Ingrese datos para la barra %0.0f: ',i);
    stringi=input(inputmessage,'s');
    material=sscanf(stringi,'%f');
    try
    Ecurr=material(1);
    Acurr=material(2);
    Icurr=material(3);
    catch
        if i==1
            error('Se necesita un primer valor para E, Area y Iz')
        else
            placeholder=0;
        end
    end

    inputmessage=sprintf('\nNodos corresp. a la barra %0.0f: ',i);
    stringi=input(inputmessage,'s');
    noditos=sscanf(stringi,'%f');
    if isempty(noditos)
        cont=false;
        break
    end
    Evec=[Evec Ecurr];
    Avec=[Avec Acurr];
    Ivec=[Ivec Icurr];
    fprintf('\n')
end

Ne=length(Evec);

