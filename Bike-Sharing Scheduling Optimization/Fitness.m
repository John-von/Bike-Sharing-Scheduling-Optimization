function FitnV = Fitness(len)
    FitnV = 1 ./ (len + 1e-10);
    FitnV(isinf(FitnV)) = 0;
    FitnV(isnan(FitnV)) = 0;
    FitnV(FitnV < 0) = 0;
end