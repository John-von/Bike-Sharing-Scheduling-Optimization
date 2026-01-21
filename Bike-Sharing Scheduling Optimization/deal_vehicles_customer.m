function [final_vehicles_customer, vehicles_used] = deal_vehicles_customer(vehicles_customer)
    mask = ~cellfun(@isempty, vehicles_customer);
    final_vehicles_customer = vehicles_customer(mask);
    vehicles_used = sum(mask);
    
    if ~isempty(final_vehicles_customer) && size(final_vehicles_customer, 2) > 1
        final_vehicles_customer = final_vehicles_customer';
    end
end