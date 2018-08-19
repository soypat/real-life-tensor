function [ ] = grafisuficiente(varargin)
%   GRAFISUFICIENTE  gráfica solicitaciones normales, cortes y momentos
%   sobre una viga de varios/un elemento(s).
%   GRAFISUFICIENTE(vectorLongitudes() ,  vectorAngulos() ,  structFuerzas{} )

% Ltotal=sum(Lv);
Lv=varargin{1};
Fv=varargin(2:end);
Ne=length(Lv);

L=[0; cumsum(Lv)];
colorN=[26, 0, 160]/256;
colorV=[133, 0, 163]/256;
colorM=[17, 147, 0]/256;

if   length(Lv)~=length(Fv) %|| length(Lv)~=length(phiv)
    error('Inputs no son de la misma cantidad de elementos.')
end

N=zeros(2*Ne,1);
V=zeros(2*Ne,1);
M=zeros(2*Ne,1);

xv=zeros(2*Ne,1);
for i=1:Ne %Generación de vectores [y1 y2] para plot
    N([i*2-1; i*2])=[-Fv{i}(1),Fv{i}(4)];
    V([i*2-1; i*2])=[Fv{i}(2),-Fv{i}(5)];
    M([i*2-1; i*2])=[-Fv{i}(3),Fv{i}(6)];
    xv([i*2-1;i*2])=[L(i) L(i+1)];
end
figure('Name','Fuerzas')

subplot(3,1,1)
    area(xv,N,'FaceColor',colorN)
    title('Fuerzas normales','fontsize',12)
    xlabel('X[in]')
    ylabel('N[lbf]')
subplot(3,1,2)
    area(xv,V,'FaceColor',colorV)
    title('Fuerzas de corte','fontsize',12)
    xlabel('X[in]')
    ylabel('V[lbf]')
subplot(3,1,3)
    area(xv,M,'FaceColor',colorM)
    title('Momentos','fontsize',12)
    xlabel('X[in]')
    ylabel('M[lbf.in]')
end

