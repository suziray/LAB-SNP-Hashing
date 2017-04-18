function CCC0s = getCCC0s(SNPs) %N samples, n SNPs

[~,n] = size(SNPs);
n = n / 2;
CCC0s = zeros(1, n*(n-1)*3);
k = 1;
for i = 1:n-1
    for j = i+1:n
        SNP1 = SNPs(:,2*i-1:2*i);
        SNP2 = SNPs(:,2*j-1:2*j);
        CCC0s(1, k:k+5) = getCCC0([SNP1 SNP2]);
        k = k + 6;
    end
end
