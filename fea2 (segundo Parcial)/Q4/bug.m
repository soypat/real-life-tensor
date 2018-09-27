syms x y real

H=0.0005;
y=H;

Tsym = ...
(204800*y)/7 - (19342813113834066795298816*y^2)/661131307601750375 + 60;

Tfunc=@(y) ...
(204800*y)/7 - (19342813113834066795298816*y^2)/661131307601750375 + 60;

symval=subs(Tsym);

funcval=Tfunc(H);

if funcval~=symval
    fprintf('Something is wrong. \n')
end


