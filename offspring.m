%%
function [new_X, new_Y] = offspring(X,Y,g,alpha)
    pop_size = numel(X);
    
    %fitness = max(1 - alpha*(g - Y(1,:))/g,0);
    fitness = (1-.001).^(g - Y(1,:));
    
    population = 1:pop_size;
    y= randsample(population,pop_size,true,fitness);
    
    new_X = X(y);
    new_Y = Y(:,y);   
end