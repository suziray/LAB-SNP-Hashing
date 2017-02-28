function SNP = randSNPGenerator(N, fa, fb)

SNP1 = rand(N,2);
SNP2 = rand(N,2);
for i = 1:N
    for j = 1:2
        if SNP1(i,j) < fa
            SNP1(i,j) = 1;
        else 
            SNP1(i,j) = 0;
        end
        if SNP2(i,j) < fb
            SNP2(i,j) = 1;
        else 
            SNP2(i,j) = 0;
        end
    end
end
SNP = [SNP1 SNP2];
end