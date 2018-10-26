%% Vigas 3D
%% Matriz global
for e=1:nelem
    long=cnod(elem(e,2),:)-cnod(elem(e,1),:);
    X=A(e)*E(e)/long;
    Y4=2*E(e)*Iz(e)/long;
    Y3=Y4*2;
    Y2=3*Y4/long;
    Y1=2*Y2/long;
    Z4=2*E(e)*Iy(e)/long;
    Z3=Z4*2;
    Z2=3*Z4/long;
    Z1=2*Z2/long;
    S=G(e)*K(e)/long;
    Ke=[X 0 0 0 0 0 -X 0 0 0 0 0
        0 Y1 0 0 0 Y2 0 -Y1 0 0 0 Y2
        0 0 Z1 0 -Z2 0 0 0 -Z1 0 -Z2 0
        0 0 0 S 0 0 0 0 0 -S 0 0
        0 0 -Z2 0 Z3 0 0 0 Z2 0 Z4 0
        0 Y2 0 0 0 Y3 0 -Y2 0 0 0 Y4
        -X 0 0 0 0 0 X 0 0 0 0 0
        0 -Y1 0 0 0 -Y2 0 Y1 0 0 0 -Y2
        0 0 -Z1 0 Z2 0 0 0 Z1 0 Z2 0
        0 0 0 -S 0 0 0 0 0 S 0 0
        0 0 -Z2 0 Z4 0 0 0 Z2 0 Z3 0
        0 Y2 0 0 0 Y4 0 -Y2 0 0 0 Y3];
    