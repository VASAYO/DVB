function PRBS = GenPRBS(N)

Reg = ones(1, 11);

PRBS = zeros(1, N);

for k = 1:N
    PRBS(k) = Reg(11);
    B = mod(sum(Reg([9, 11])), 2);
    Reg(2:11) = Reg(1:10);
    Reg(1) = B;
end