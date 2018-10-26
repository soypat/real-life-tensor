function Disp = Solve_bar(Ab,At,L,P,nel,A_type)
% Discretizacion
nnod = nel + 1;                  %Numero de Nodos
Element = zeros(2,nel);
areaNod = linspace(Ab,At,nnod);

% Area mayor
%     A = areaNod(1:end-1);
%
%     %Area menor
%     A = areaNod(2:end);
%
if strcmp(A_type,'P')
    %Area promedio
    A = filter(ones(1,2)/2,1,areaNod);
    A(1) = [];
elseif strcmp(A_type,'m')
    %Area menor
    A = areaNod(2:end);
elseif strcmp(A_type,'M')
    %Area mayor
    A = areaNod(1:end-1);
end

Node = (0:L/nel:L)';             %Coordenadas Nodales
Element(1,:) = 1:nnod - 1;
Element(2,:) = 2:nnod;
Element = Element';              %Vector Conexion Elementos
%
Bc = zeros(nnod,1);
Bc(1) = 1;                       %Vector Condiciones de Borde
R = zeros(nnod,1);
R(nnod) = P;                 %Vector Cargas

% Propiedades del Material
E = 210000;

% Armado Matriz de Rigidez
K = zeros(nnod,nnod);
for e = 1:nel;
    Node_el = Element(e,1:2);
    k = A(e)*E/abs(Node(Node_el(2)) - Node(Node_el(1)));
    Ke = [ k -k
        -k  k ];
    K(Node_el,Node_el) = K(Node_el,Node_el) + Ke;
end

% Reduccion Matriz
e = find(Bc);
R(e) = [];
K(e,:) = [];
K(:,e) = [];

% Solver
D = K\R;

% Reconstruccion
e = find(Bc - 1);
De = zeros(nnod,1);
De(e) = De(e) + D;

% Tensiones
Se = zeros(1,nel);
for e = 1:nel;
    Node_el = Element(e,:);
    B = [-1 1]/abs(Node(Node_el(2)) - Node(Node_el(1)));
    Se(e) = dot(B,De(Node_el))*E;
end
Disp = De(nnod);

% Salida de datos

%     disp(['Número de elementos: ' num2str(nel)])
%
%     disp('Configuracion inicial')
%     disp(Node)
%
%     disp('Desplazamientos')
%     disp(De)
%
%     disp('Configuración deformada')
%     disp(De + Node)
%
%     disp('Tensiones normales')
%     disp(Se')
