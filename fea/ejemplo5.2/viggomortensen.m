nodos=12*[-5,0;-5,10;5,10;5,0];%in
E=30e6;%psi
A=10;%in^2
I=200;%in^4

viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];
vigota=@(E,A,I,L) [A*E/L 0 0 -A*E/L 0 0;
    0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
    0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
    -A*E/L 0 0 A*E/L 0 0;
    0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -E*I/L^2;
    0 6*E*I/L^2 2*E*I/L -6*E*I/L^2 4*E*I/L];

Tv=@(phi) [cosd(phi) sind(phi) 0 0 0 0;
    -sind(phi) cosd(phi) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cosd(phi) sind(phi) 0;
    0 0 0 -sind(phi) cosd(phi) 0;
    0 0 0 0 0 1];
