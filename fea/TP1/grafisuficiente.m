function [ ] = grafisuficiente(varargin)
%   GRAFISUFICIENTE  gr�fica solicitaciones normales, cortes y momentos
%   sobre una viga de varios/un elemento(s).
%   GRAFISUFICIENTE(vectorLongitudes() ,  vectorAngulos() ,  structFuerzas{} )

% Ltotal=sum(Lv);
Lv=varargin{1};
Fcell=varargin(3:end);
phiv=varargin{2};
Ne=length(Lv);

L=[0; cumsum(Lv)];
colorN=[26, 0, 160]/256;
colorV=[133, 0, 163]/256;
colorM=[17, 147, 0]/256;



if   length(Lv)~=length(Fcell) || length(Lv)~=length(phiv)
    error('Inputs no son de la misma cantidad de elementos.')
end

N=zeros(2*Ne,1);
V=zeros(2*Ne,1);
M=zeros(2*Ne,1);

xv=zeros(2*Ne,1);
for i=1:Ne %Generaci�n de vectores [y1 y2] para plot
    %ROTACION DE VECTOR
    s=sind(-phiv(i));
    c=cosd(-phiv(i));
    Fv=Fcell{i};
    Fv([1 2 4 5])=[Fv(1)*c-Fv(2)*s;c*Fv(2)+s*Fv(1);Fv(4)*c-Fv(5)*s;c*Fv(5)+s*Fv(4)];
    %creaci�n vectores
    N([i*2-1; i*2])=[-Fv(1),Fv(4)];
    V([i*2-1; i*2])=[Fv(2),-Fv(5)];
    M([i*2-1; i*2])=[-Fv(3),Fv(6)];
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
