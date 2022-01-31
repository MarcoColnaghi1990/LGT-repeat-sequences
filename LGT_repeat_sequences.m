%%% Simulation Parameters
g = 100;                                    %%% Genome size
N = 2500;                                   %%% Population size
rho = 0.1;                                  %%% Rep Seq density
t_max = 5000;                               %%% N generations
alpha = 1;                                  %%% Slope of fitness curve
n_runs = 100;                               %%% Number of replicates
u = 3e-05;                                  %%% mutation rate
lambda = 0.1;                               %%% LGT rate
L = 10;                                     %%% Recombination length

%%% Initialise variables
meanGeneLoss = zeros(1,n_runs);
meanMutLoad = zeros(1,n_runs);
fixDel = zeros(1,n_runs);
fixMut = zeros(1,n_runs);


%% Mutation accumulation in the presence of repeats


    for n = 1:n_runs
        [meanGeneLoss(n), meanMutLoad(n), fixDel(n), fixMut(n)] = ...
            EvolutionWithRepSeqs_Mutations_mu(g,N,L,lambda, rho, t_max,alpha,2,u);
    end

    
average_gene_loss_deletion = mean(meanGeneLoss/t_max);
average_gene_loss_mutation = mean(meanMutLoad/t_max);
