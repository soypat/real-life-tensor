X = E*A/Le(i);
            Y1 = 12*E*Iz/Le(i)^3;
            Y2 = 6*E*Iz/Le(i)^2;
            Y3 = 4*E*Iz/Le(i)^1;
            Y4 = 2*E*Iz/Le(i)^1;
            Kdiag = diag([X Y1 Y3 X Y1 Y3]);
            Kp = [0     0       0       -X      0       0
                  0     0       Y2      0       -Y1     Y2
                  0     0       0       0       -Y2     Y4
                  0     0       0       0       0       0
                  0     0       0       0       0       -Y2
                  0     0       0       0       0       0];
            Kel = Kp + Kp' + Kdiag;