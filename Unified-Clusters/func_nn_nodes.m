%% Neural Gas Nearest Neighbor Fidner
function nn_dist = func_nn_nodes(nodes, w)

for i = 1: nodes
%     min_dist = 10^23; 
    distanc = [];
    for j = 1:nodes
        ds2 = 0;
        for x = 1:numel(w(1,:))
            ds2 = ds2 + (w(i,x) - w(j,x))^2;
        end
        distanc(j) = ds2^0.5;

    end
    nn_dist(i,1) = min(distanc(distanc >0));
end

end