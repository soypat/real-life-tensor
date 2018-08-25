function [ ] = grafisuficiente2(varargin)
%   GRAFISUFICIENTE  gráfica solicitaciones normales, cortes y momentos
%   sobre una viga de varios/un elemento(s).
%   GRAFISUFICIENTE(vectorLongitudes() ,  structFuerzas{} )

% Ltotal=sum(Lv);

Lv=varargin{1};
Fcell=varargin(2:end);
Ne=length(Lv);

L=[0; cumsum(Lv)];
colorN=[26, 0, 160]/256;
colorV=[133, 0, 163]/256;
colorM=[17, 147, 0]/256;



if   length(Lv)~=length(Fcell)% || length(Lv)~=length(phiv)
    error('Inputs no son de la misma cantidad de elementos.')
end

N=zeros(2*Ne,1);
V=zeros(2*Ne,1);
M=zeros(2*Ne,1);

xv=zeros(2*Ne,1);
for i=1:Ne %Generación de vectores [y1 y2] para plot
    %ROTACION DE VECTOR
    Fv=Fcell{i};
%     Fv([1 2 4 5])=[Fv(1)*c-Fv(2)*s;c*Fv(2)+s*Fv(1);Fv(4)*c-Fv(5)*s;c*Fv(5)+s*Fv(4)];
    %creación vectores
    for j=1:length(Fv)
        if abs(Fv(j))<1e-9
            Fv(j)=0;
        end
    end
    if length(Fv)==2
        N([i*2-1; i*2])=[-Fv(1),Fv(2)];
        continue;
    end
    N([i*2-1; i*2])=[-Fv(1),Fv(4)];
    V([i*2-1; i*2])=[Fv(2),-Fv(5)];
    M([i*2-1; i*2])=[-Fv(3),Fv(6)];
    xv([i*2-1;i*2])=[L(i) L(i+1)];
end

figure('Name','Fuerzas2')

subplot(3,1,1)
    area(xv,N,'FaceColor',colorN)
    title('Fuerzas normales','fontsize',12)
    xlabel('X[in]')
    ylabel('N[lbf]')
    xlim([0 L(end)])
subplot(3,1,2)
    area(xv,V,'FaceColor',colorV)
    title('Fuerzas de corte','fontsize',12)
    xlabel('X[in]')
    ylabel('V[lbf]')
    xlim([0 L(end)])
subplot(3,1,3)
    area(xv,M,'FaceColor',colorM)
    title('Momentos','fontsize',12)
    xlabel('X[in]')
    ylabel('M[lbf.in]')
    xlim([0 L(end)])
end

