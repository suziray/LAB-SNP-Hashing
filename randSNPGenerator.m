function SNP = randSNPGenerator(N, fa)

SNP = rand(N,2);
for i = 1:N
    for j = 1:2
        if SNP(i,j) < fa
            SNP(i,j) = 1;
        else 
            SNP(i,j) = 0;
        end
    end
end
end