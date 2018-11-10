if k*(smin+x0)>m*abs(amin)
    fprintf('\n1...')
    cond=1;
else
    fprintf('\n\nCondicion 1 fallida\nk*(smin+x0)>m*abs(amin)\n\n')
    cond=0;
end
if min(Fc)>0 && min(Fc)>0.09*max(Fc)
    cond=cond+1;
    fprintf('2...')
else
    fprintf('\nCondicion 2 fallida. Fc<0 OR Fcmin>0.1*Fcmax!\n\n')
end

if w/w_resonancia<sqrt(2)/2
    cond=cond+1;
    fprintf('3...\n%0.0f/3 Condiciones cumplidas',cond)
else
    fprintf('\n\nCondicion resonancia fallida.\n\n%0.0f/3 condiciones cumplidas\n',cond);
end