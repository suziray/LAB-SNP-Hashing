N = 1000;
fA = 0.7;
fB = 0.9;
q = 1.5;

SNP1 = rand(N,2);
SNP2 = rand(N,2);
for i = 1:N
    for j = 1:2
        if SNP1(i,j) > fA
            SNP1(i,j) = 1;
        else
            SNP1(i,j) = 0;
        end
        if SNP2(i,j) > fB
            SNP2(i,j) = 1;
        else
            SNP2(i,j) = 0;
        end
    end
end

R = zeros(2,2);
for i = 1:N
    if SNP1(i,1)*SNP1(i,2) == 0
        if SNP2(i,1)*SNP2(i,2) == 0
            R(1,1) = R(1,1) + 1;
        else
            R(1,2) = R(1,2) + 1;
        end
    else
         if SNP2(i,1)*SNP2(i,2) == 0
            R(2,1) = R(2,1) + 1;
        else
            R(2,2) = R(2,2) + 1;
         end
    end
end
R = R/N;

F = zeros(2,2);
for i = 1:N
    for j = 1:2
        if SNP1(i,j) == 0
            F(1,1) = F(1,1) + 1;
        else
            F(1,2) = F(1,2) + 1;
        end
        if SNP2(i,j) == 0
            F(2,1) = F(2,1) + 1;
        else
            F(2,2) = F(2,2) + 1;
        end
    end
end
F = F / 2 / N;
F = 1 - F / q;

CCC = zeros(2,2);
for i = 1:2
    for j = 1:2
        CCC(i,j) = 9 / 2 * R(i,j) * F(1,i) * F(2,j);
    end
end
