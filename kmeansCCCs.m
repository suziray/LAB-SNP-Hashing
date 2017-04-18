function idx = kmeansCCCs(SNPs, k) %N samples, n SNPs

SNPs = transpose(SNPs);
idx2 = kmeans(SNPs, k, 'Distance', 'correlation');
[n,~] = size(SNPs);
n = n / 2;
idx = zeros(1, n);

for i = 1:n
    id1 = idx2(2*i - 1, 1);
    id2 = idx2(2*i, 1);
    if id1 == id2
        idx(1, i) = id1;
    else
        sum1 = 0; n1 = 0;
        v1 = SNPs(2*i - 1, :);
        v2 = SNPs(2*i, :);
        sum2 = 0; n2 = 0;
        for j = 1:2*n
            if idx2(j, 1) == id1 && j ~= 2*i - 1
                sum1 = sum1 + corrcoef(v1,SNPs(j, :));
                n1 = n1 + 1;
            elseif idx2(j, 1) == id2 && j ~= 2*i
                sum2 = sum2 + corrcoef(v2,SNPs(j, :));
                n2 = n2 + 1;
            end
        end
        if sum1/n1 > sum2/n2
            idx(1, i) = id1;
        else
            idx(1, i) = id2;
        end
    end
end