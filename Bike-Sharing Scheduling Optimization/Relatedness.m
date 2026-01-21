function Rij = Relatedness(i, j, dist, vehicles_customer)
    NV = size(vehicles_customer, 1);
    d = dist(i+1, j+1);
    [md, ~] = max((dist(i+1, 2:end)));
    c = d / md;
    V = 1;
    for k = 1:NV
        route = vehicles_customer{k};
        findi = find(route == i, 1, 'first');
        findj = find(route == j, 1, 'first');
        if ~isempty(findi) && ~isempty(findj)
            V = 0;
        end
    end
    Rij = 1 / (c + V);
end