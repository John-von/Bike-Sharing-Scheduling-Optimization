function ObjV = calObj_hybrid(Chrom, cusnum, cap, demands, a, b, L, s, dist, alpha, beta, user_cost, gamma, operator_idx)
    NIND = size(Chrom, 1);
    ObjV = zeros(NIND, 1);
    for i = 1:NIND
        [VC, ~, ~, ~, ~] = decode_hybrid(Chrom(i,:), cusnum, cap, demands, a, b, L, s, dist, operator_idx);
        costF = costFuction(VC, a, b, s, L, dist, demands, cap, alpha, beta);
        ObjV(i) = costF + gamma * user_cost;
    end
end
