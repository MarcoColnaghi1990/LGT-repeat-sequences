%%
function [X,Y] = Recombination_With_Repeats(X,X_old,L,N,n_genes,sigma,Y)
i = randi(N);
g = numel(X_old{i});
r = randi(g-L+1);                   %%% Select random locus
eDNA_seq = X_old{i}(r:r+L-1);      %%% eDNA sequence from previous generation


seq_1 = eDNA_seq(1);                %%% Flanking sequence 1
seq_2 = eDNA_seq(end);              %%%    "         "    2

seq_1_position=find(abs(X)==abs(seq_1));   %%% Find homologs of seq_1

if numel(seq_1_position)>1
    rec_start = randsample(seq_1_position,1);       %%% beginning of recombination
else
    rec_start = seq_1_position;
end

new_list = [X(rec_start:end),X(1:rec_start-1)];     %%% circular chromosome


h =  find(abs(new_list) == abs(seq_2));            %%% Loci with homologous to seq_2


if numel(h)>0
    weight = 1./exp((h-L).^2/(2*sigma*L));    %%% Recombination probability
    if numel(h)>1 && sum(weight)>0
        
        y = randsample(1:numel(h),1,true,weight);
    else
        y=1;
    end
    if rand()<weight(y)
        deleted_seq = new_list(1:h(y));
        if numel(eDNA_seq) ~= h(y)
            new_list = [eDNA_seq, new_list(h(y)+1:end)];
            Y(1) = intersect2(X,Y(1),eDNA_seq,deleted_seq);
            X = new_list;
            %Y(1) = numel(intersect(1:n_genes,X));
            Y(2) = numel(X);
        else
            if sum(eDNA_seq ~= new_list(1:h(y)))>0
                Y(1) = intersect2(X,Y(1),eDNA_seq,deleted_seq);
                new_list = [eDNA_seq, new_list(h(y)+1:end)];
                X = new_list;
                %Y(1) = numel(intersect(1:n_genes,X));
                Y(2) = numel(X);
                %disp('It happened now!')
            end
        end
    end
end

end


%%
function Y = intersect2(X,Y,list1,list2)

%%% eliminate repeats
list2(list2==0)=[];
list1(list1==0)=[];

for i = 1:numel(list1)
    a = list1(i);
    if a > 0                    %%% wild-type in new sequence
        if sum(list2==a)>0      %%% if it's in deleted seq -> no change
        else                    %%% if it's not, is it in the genome?
            if sum(X==a) == 0
                Y = Y + 1;      %%% gene gain
            end
        end
    else                        %%% mutant in the new sequence
        number_of_copies = sum(list2==abs(a));    %%% of wild-type in the deleted sequence
        if number_of_copies>0   %%% wild-type in the old sequence
            if sum(X==abs(a)) == number_of_copies       %%% not elsewhere in the genome
                Y = Y - 1;      %%% gene loss
            end
        end
    end
end

for i = 1:numel(list2)
    b = list2(i);
    if b > 0        %%% wild-type in deleted sequence
        if sum(abs(list1)==b) == 0       %%% not present/mutated in new sequence
            Y = Y - 1;          %%% gene loss
        end
    end
end

end





%%% For each gene in list 1 there are
%%% 6 possibilities:

%%% the wild-type gene is in the new seq, wild-type in the del seq
%%%NO CHANGES

%%% the wild-type gene is in the new seq, mutant in the del seq
%%%GENE GAIN (if the gene is not elsewere in the genome)
%%%NO CHANGES (otherwise)

%%% the wild-type gene is in the new seq, not in the del seq
%%%GENE GAIN (if the gene is not elsewere in the genome)
%%%NO CHANGES (otherwise)

%%% the mutant gene is in the new seq, wild-type in the del seq
%%%GENE LOSS (if the gene is not elsewere in the genome)
%%%NO CHANGES (otherwise)

%%% the mutant gene is in the new seq, mutant in the del seq
%%%NO CHANGES

%%% the mutant gene is in the new seq, not in the del seq
%%%NO CHANGES



%%% For each gene in list 2 there are
%%% 2 possibilities:

%%% the wild-type gene is in the del seq, not in the new seq
%%%GENE LOSS (if the gene is not elsewere in the genome)
%%%NO CHANGES (otherwise)

%%% the mutant gene is in the del seq, not in the new seq
%%%NO CHANGES
