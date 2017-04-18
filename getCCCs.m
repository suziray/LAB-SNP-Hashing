function CCCs = getCCCs(SNPs) %N samples, n SNPs

[~,n] = size(SNPs);
n = n / 2;
CCCs = zeros(1, n*(n-1)/2);
k = 1;
for i = 1:n-1
    for j = i+1:n
        SNP1 = SNPs(:,2*i-1:2*i);
        SNP2 = SNPs(:,2*j-1:2*j);
        CCCs(1, k) = getCCC([SNP1 SNP2]);
        k = k + 1;
    end
end
