function y = Pol(x,coef)

    tamanho = numel(coef);
    
    resp = 0;
    for i = 1:tamanho
        
        n = i-1;
        termo = coef(i);
        resp = resp+termo*x^n; 
    end
    y = resp;