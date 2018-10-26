%% Mallador axisimetrico
clear; close all; clc

%% Geometría
a = 4;
b = 8;
d = 2;

%% Tipo de elemento y tamaño
eleT = 'Q8'; 
SZ = 2;
SR = 4;
nel = SR*SZ;

%% Nodos
switch eleT
    case 'Q4'
        pZ = d/SZ;
        pR = (b-a)/SR;
        [locR,locZ] = meshgrid(a:pR:b,0:pZ:d);
        nodes = [reshape(locR',[],1) reshape(locZ',[],1)];
    case 'Q8'
        pZ = d/SZ/2;
        pR = (b-a)/SR/2;
        [locR,locZ] = meshgrid(a:pR:b,0:pZ:d);
        nodes = [reshape(locR',[],1) reshape(locZ',[],1)];
        keep = true(size(nodes));
        for iz = 1:SZ
            for ir = 1:SR
                s = (1+2*SR)*(2*iz-1) + ir*2;
                keep(s,:) = false;
            end
        end
        nodes = [nodes(keep(:,1),1) nodes(keep(:,2),2)];
    case 'Q9'
        pZ = d/SZ/2;
        pR = (b-a)/SR/2;
        [locR,locZ] = meshgrid(a:pR:b,0:pZ:d);
        nodes = [reshape(locR',[],1) reshape(locZ',[],1)];
end
%% Elementos
switch eleT
    case 'Q4'
        elements = [];
        base = (1:1:SR)';
        bloq = [base base+1 base+SR+2 base+SR+1];
        for iz = 1:SZ
            elements = [elements; bloq+(SR+1)*(iz-1)];
        end
    case 'Q8'
        elements = [];
        baseA = (1:2:2*SR-1)';
        baseB = (1:1:SR)';
        bloq = [baseA         baseA+2      baseA+3*SR+4      baseA+3*SR+2 ... 
                baseA+1    baseB+2*SR+2    baseA+3*SR+3      baseB+2*SR+1];
        for iz = 1:SZ
            elements = [elements; bloq+(3*SR+2)*(iz-1)];
        end
end

% meshplot(elements,nodes,'b')
nodename = strcat(['Nodos',num2str(SZ),'x',num2str(SR),num2str(eleT),'.mat']);
elename = strcat(['Elementos',num2str(SZ),'x',num2str(SR),num2str(eleT),'.mat']);
save(nodename,'nodes');
save(elename,'elements');