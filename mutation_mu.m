%%
function [X,Y] = mutation_mu(X,Y,N,u)

random_list = rand(1,N);
U = u*Y(2,:);
list_1 = find(random_list<U);

for k=1:numel(list_1)
    i=list_1(k);
    genes_positions=1:numel(X{i});%genes_positions = find(X{i}~=0);
    y = randsample(genes_positions,1);
    %for each mutated gene, calculate the number of copies:
    if X{i}(y)>0
        if(numel(X{i}==X{i}(y)==1))
            Y(1,i)=Y(1,i)-1;
        end
    end
    X{i}(y) = -abs(X{i}(y));
    
end


end