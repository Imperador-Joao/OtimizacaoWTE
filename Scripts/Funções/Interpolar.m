function y0 = Interpolar(x,y,x0)

    %Verificar se os tamanhos são iguais

    if numel(x) ~= numel(y)
        error("Tamanhos incompatíveis")
    end
    
    %Se y for um vetor linha, ele transformado em um vetor coluna
    
    if size(y,2) > size(y,1)
        y = y';
    end
    
    % Gerar um polinômio de grau 2 em torno de (x0,y0)
    
    x_prox = min(abs(x-x0))+x0;                                            %Valor de x mais próximo de x0 para mais
    if ~Pertence_a(x_prox,x)
        x_prox = -min(abs(x-x0))+x0;                                       %Valor de x mais próximo de x0 para menos           
    end
    
    pos = find(x==x_prox);
    y_prox = y(pos);
    
    if pos > 1
    
        x_antes = x(pos-1);
        y_antes = y(pos-1);

        x_depois = x(pos+1);
        y_depois = y(pos+1);


        x = [x_antes,x_prox,x_depois];                                     %Substitui x e y por uma versão reduzida com 3 pontos.
        y = [y_antes;y_prox;y_depois];
        
    else
        
        x_depois = x(pos+1);
        y_depois = y(pos+1);
        
        x_depois_depois = x(pos+2);
        y_depois_depois = y(pos+2);
        
        x = [x_prox;x_depois;x_depois_depois];
        y = [y_prox;y_depois;y_depois_depois];
    end
    
    %Resolver o sistema M*x = y
     
    tamanho = numel(x);
    Matriz_sistema = [];
    
    for i = 1:tamanho
        
        linha = [];
        xn = x(i);
        
        for j = 1:tamanho
            
           n = j-1;
           linha = [linha xn^n];
        
        end
        
        Matriz_sistema = [Matriz_sistema;linha];
        
    end
    
    polinomio = Matriz_sistema\y;
    y0 = Pol(x0, polinomio);
    
    
    