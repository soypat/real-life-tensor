%% L'Arc (como el boliche)
%Lo voy a hacer tomando cosas del ejercicio 0, osea, "de cero"
clear dN N dNaux %Esto para evitar los bugs comunmente ocasionados por el symengine de matl�
%% Defino el elemento a usar (Q8) y hallo las funciones de forma
%Q8 t�pico con numeraci�n del Cook y de la catedra.
syms ksi eta real
% X = [1 ksi eta ksi^2 eta^2 ksi*eta ksi^2*eta eta^2*ksi];
X = [1 ksi eta ksi^2 eta^2 ksi*eta ksi^2*eta eta^2*ksi];
A = [1 -1 -1  1  1  1 -1 -1 %Nodo 1
     1  1 -1  1  1 -1 -1  1 % 2
     1  1  1  1  1  1  1  1 % 3
     1 -1  1  1  1 -1  1 -1  %4 
     1  0 -1  0  1  0  0  0 %5
     1  1  0  1  0  0  0  0 %6
     1  0  1  0  1  0  0  0 %7
     1 -1  0  1  0  0  0  0]; %Nodo 8  (a manopla esto chicos, deberia ser como respirar)
shapefuns = X/A;
%Igual que el ejercicio 0
N(1,1:2:2*length(shapefuns))=shapefuns;
N(2,2:2:2*length(shapefuns))=shapefuns; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). Despues veo si me sirven
dN(1,1:2:2*length(shapefuns))=diff(shapefuns,ksi);
dN(2,2:2:2*length(shapefuns))=diff(shapefuns,eta);
dNaux=[diff(shapefuns,ksi);diff(shapefuns,eta)]; %Para calcular jacobiano
%% Carga Malla ADINA (ya mallado en Q8, duh)
aux=load('nod.txt');
nodos=aux(:,2:3);
aux=load('elem_arc.txt');
elem=aux(:,2:9);%porque es Q8
%% Datos Material
E=200e9; %Trabajar unidades consistentes tipo metros,Pascales,Newtons
nu=0.3;
t=1;

%% DOF Define, the Dofinitions ;)
[Nnod, Ndofpornodo]=size(nodos); %N mayuscula=N�mero
[Nelem, nodporelem]=size(elem);
doftot=Nnod*Ndofpornodo;
dof = reshape(1:doftot,Ndofpornodo,Nnod)';
%% Gauss orden 3
npg=9;
a=sqrt(.6);
upg=[-a -a
     -a  0
      -a a
      0  -a
      0  0
      0  a
      a  -a
      a  0
      a  a];
  wpg=5/9*zeros(npg,1);
  wpg(5)=8/9;
%% Rigidity Matrix. 
% ukno i lik em rigid
Cstress=E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 .5-nu]; %Constitutive 2D stressful
K=zeros(doftot);
tic
for e=1:Nelem
    Ke=zeros(Ndofpornodo*nodporelem);
    index=elem(e,:);
    elenod=nodos(index,:);
    for ipg=1:npg
        ksi=upg(ipg,1); eta=upg(ipg,2);
        
        J=double(subs(dNaux))*elenod;
        
        dNxy=J\double(subs(dN)); %El double para que tarde menos, hace que se saltee la algebra (resultado exacto)
        
        B=zeros(size(Cstress,2),Ndofpornodo*nodporelem);
        B(1:2,:)=dNxy;
        B(3,1:2:end)=dNxy(2,2:2:end)
        B(3,2:2:end)=dNxy(2,1:2:end)
        
        B2=zeros(size(Cstress,2),size(Ke,1)); %Segunda Forma
        B2(1:2,:)=dNxy;
        B2(3,2:end)=dNxy(1,1:end-1);
        B2(3,1:end-1)=B(3,1:end-1)+dNxy(2,2:end);
    end
end
toc