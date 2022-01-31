%%
function [X1, X2, X3, X4] = EvolutionWithRepSeqs_Mutations_mu (g, N, L, lambda, rho, t_max, alpha,sigma,U)
%% Initialisation
R_vector=round(rho*g)*ones(1,N);
X = cell(1,N);
nrepeats = round(rho*g);
Y = g*ones(2,N);
Y(2,:)=Y(2,:)+nrepeats;
y = sort(randsample(1:g+nrepeats,g));

for i = 1:N
    X{i}=zeros(1,g+R_vector(i));
    X{i}(y)=1:g;
end

num_genes=zeros(1,N);information_degradation=zeros(1,N);
num_mutations=zeros(1,N);

%% Evolutionary dynamics
for time = 1:t_max
    X_old=X;
    [X,Y] = offspring(X,Y,g,alpha);
    [X,Y] = mutation_mu(X,Y,N,U);
    random_list = rand(1,N);
    for i = 1:N
        if random_list(i) < lambda
            [X{i},Y(:,i)] = Recombination_With_Repeats(X{i},X_old,L,N,g,2,Y(:,i));
        end
    end
    
end

for i = 1:N
    list = X{i}(X{i}~=0);
    genes = numel(intersect(abs(list),1:g));
    wildtypes = numel(intersect(list,1:g));
    num_genes(i)=genes;  %%%total number of genes (mutated + wild-type)
    information_degradation(i) = g - num_genes(i); %%%deletions
    num_mutations(i) = genes-wildtypes; %%%
end

X1 = mean(information_degradation);     %%% Number of different wild-type genes
X2 = mean(num_mutations);               %%% Average mut load
X3 = min(information_degradation);      %%% Least Loaded Class (Deletions)
X4 = min(num_mutations);                %%% Least Loaded Class (Mutations)
end