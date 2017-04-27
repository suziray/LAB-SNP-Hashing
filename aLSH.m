function [h0, h1] = aLSH(SNPs)

[d, n] = size(SNPs);
m = 3;
U = 0.83;
r = 2.5;
N = 100;

nSNPs = scaling(SNPs, U, n);
qs = zeros(d + 2 * m, n);
for i = 1:n
    qs(:,i) = getQ(nSNPs(:,i), m ,d);
end
ps = zeros(d + 2 * m, n);
for i = 1:n
    ps(:,i) = getP(nSNPs(:,i), m, d);
end

for i = 1:N
    a = normrnd(0, 1, [1 d+2*m]);
    b = r * rand;
    h0 = zeros(1, n);
    h1 = zeros(1, n);
    for j = 1:n
        h0(1,j) = hash(qs(:,j),a,b,r);
        h1(1,j) = hash(ps(:,j),a,b,r);
    end
end

end

function SNP1s = scaling(SNPs, U, n)
M = norm(SNPs(:,1));
for i = 2:n
    M = max(M, norm(SNPs(:,i)));
end
SNP1s = SNPs * U / M;
end

function q1 = getQ(q, m, d)
q1 = ones(d + 2 * m, 1) / 2;
q1(1:d,1) = q;
q1(d+m+1,1) = sqrt(norm(q));
for i = 2:m
    q1(d+m+i,1) = sqrt(q1(d+m+i-1,1));
end
end

function p1 = getP(x, m, d)
p1 = ones(d + 2 * m, 1) / 2;
p1(1:d,1) = x;
p1(d+1,1) = sqrt(norm(x));
for i = 2:m
    p1(d+i,1) = sqrt(p1(d+i-1,1));
end
end

function h = hash(x, a, b, r)
h = floor((a*x+b)/r);
end
