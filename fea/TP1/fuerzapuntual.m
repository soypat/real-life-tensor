function [R] = fuerzapuntual(R,minodo,fx,fy,mz)
vec=[3*minodo-2 3*minodo-1 3*minodo];
R(vec)=[fx fy mz];
end

