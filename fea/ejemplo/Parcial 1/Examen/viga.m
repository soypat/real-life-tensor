%% Funcion viga
function [Kelemento] = viga(A,Iv,E,cnod,elem,e)
    v = cnod(elem(e,2),:) - cnod(elem(e,1),:);   % Vector del elemento
    long = norm(v);                              % Norma del elemento
    vd = v/long;                                 % Vector director
        X  = E*A/long;
        Y1 = 12*E*Iv/long^3;
        Y2 = 6*E*Iv/long^2;
        Y3 = 4*E*Iv/long;
        Y4 = 2*E*Iv/long;
        Kviga =[X      0       0      -X       0        0
                0     Y1      Y2       0     -Y1       Y2
                0     Y2      Y3       0     -Y2       Y4
               -X      0       0       X       0        0
                0    -Y1     -Y2       0      Y1      -Y2
               0     Y2      Y4       0     -Y2       Y3];
        lambda = [vd 0; -vd(2) vd(1) 0; 0 0 1];
        T = blkdiag(lambda,lambda);
        Kelemento = T'*Kviga*T;
end