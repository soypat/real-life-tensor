function [T] = Tbu(phi)
%TBU matriz transf. de barra
T=[cosd(phi) 0;
    sind(phi) 0;
    0 cosd(phi);
    0 sind(phi)];
end

