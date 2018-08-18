function GrafitodoE(varargin)
%GRAFITODO como lo dice la funcion, grafica todo
%ingresar primero las longitudes de los elementos y despues las fuerzas
%locales

n=nargin;
Lv=[0 varargin{1:n/2}];
L=cumsum(Lv);
fl={varargin{n/2+1:n}};

figure('Name','Esfuerzos')
if length(fl{1})==4%solo calcula corte y flexion
    %ahora grafico el corte
    for i=1:length(L)-1
        l=length(L(i):L(i+1));
        V=linspace(fl{i}(1),-fl{i}(3),l);
        M=linspace(-fl{i}(2),fl{i}(4),l);
        
        subplot(2,1,1)
        hold on
        plot(0:L(end),linspace(0,0,length(0:L(end))),'k')%grafico el cero
        plot(linspace(L(i),L(i),l),linspace(0,fl{i}(1),l),'r')%abro
        plot(L(i):L(i+1),V,'r')%grafico
        plot(linspace(L(i+1),L(i+1),l),linspace(-fl{i}(3),0,l),'r')%cierro
        title('Tensiones de corte','fontsize',12) 
        xlabel('X[in]')
        ylabel('Tc[psi]')
        %grafico el momento
        subplot(2,1,2)
        hold on
        plot(0:L(end),linspace(0,0,length(0:L(end))),'k')%grafico el cero
        plot(linspace(L(i),L(i),l),linspace(0,-fl{i}(2),l),'b')%abro
        plot(L(i):L(i+1),M,'b')%grafico
        plot(linspace(L(i+1),L(i+1),l),linspace(0,fl{i}(4),l),'b')%cierro
        title('Tensiones flexoras','fontsize',12)
        xlabel('X[in]')
        ylabel('Sf[psi]')
    end
elseif length(fl{1})==6 %calcula normales tambien, quiere decir que esta rotado
    
    for i=1:length(L)-1
        
        
        l=length(L(i):L(i+1));
        N=linspace(-fl{i}(1),fl{i}(4),l);
        V=linspace(fl{i}(2),-fl{i}(5),l);
        M=linspace(-fl{i}(3),fl{i}(6),l);
        %ahora grafico las normales
        
        subplot(3,1,1)
        hold on
        plot(0:L(end),linspace(0,0,length(0:L(end))),'k')%grafico el cero
        plot(linspace(L(i),L(i),l),linspace(0,-fl{i}(1)),'r')%abro
        plot(L(i):L(i+1),N,'r')%grafico
        plot(linspace(L(i+1),L(i+1),l),linspace(-fl{i}(4),0,l),'r')%cierro
        title('Tensiones normales','fontsize',12)
        xlabel('X[in]')
        ylabel('Sa[psi]')
        
        %ahora grafico el corte
        subplot(3,1,2)
        plot(0:L(end),linspace(0,0,length(0:L(end))),'k')%grafico el cero
        plot(linspace(L(i),L(i),l),linspace(0,fl{i}(2),l),'r')%abro
        plot(L(i):L(i+1),V,'r')%grafico
        plot(linspace(L(i+1),L(i+1),l),linspace(-fl{i}(5),0,l),'r')%cierro
        title('Tensiones de corte','fontsize',12)
        xlabel('X[in]')
        ylabel('Tc[lbf]')
        %grafico el momento
        subplot(3,1,3)
        hold on
        plot(0:L(end),linspace(0,0,length(0:L(end))),'k')%grafico el cero
        plot(linspace(L(i),L(i),l),linspace(0,-fl{i}(3),l),'b')%abro
        plot(L(i):L(i+1),M,'b')%grafico
        plot(linspace(L(i+1),L(i+1),l),linspace(0,fl{i}(6),l),'b')%cierro
        title('Tensiones flexoras','fontsize',12)
        xlabel('X[in]')
        ylabel('Sf[psi]')
    end
end

end








