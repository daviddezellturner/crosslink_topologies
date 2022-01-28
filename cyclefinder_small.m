
function [min_cycle, min_cost] = cyclefinder_small(graph)

    num_nodes = size(graph,1);

    cycles = cell(0);
    % only allow one revisit for now?
    combo_set = [1:num_nodes];

    base = 2:num_nodes;
    
    for extra_length = 0:num_nodes;
        nodeset_new = combntns(combo_set,extra_length);
        size_nn = size(nodeset_new);
        num_nn = size_nn(1);

        for j=1:num_nn
            extra_nodes = nodeset_new(j,:);
            all_nodes = [base extra_nodes];
            p = perms(all_nodes);
            size_p = size(p);
            num_p = size_p(1);
            for k=1:num_p
                cycle = [1 p(k,:) 1];
                cycles{end+1} = cycle;
            end
        end

    end

   size_cycles = size(cycles);
   num_cycles = size_cycles(2);
   
   min_cost = -1;
   min_cycle = [];
   for i = 1:num_cycles
       cycle = cycles(i);
       cycle = cycle{1};
       cost = -1;
       bottleneck = -1;
       for j = 1:length(cycle)-1
           if graph(cycle(j),cycle(j+1)) ~= 0
               edge_cost = graph(cycle(j),cycle(j+1));
               if edge_cost > bottleneck
                   bottleneck = edge_cost;
               end
           else
               bottleneck = -1;
               cost = -1;
               break
           end
       end
       if bottleneck ~= -1
           cost = sum(1:length(cycle)-1)*bottleneck;
       end
       if cost ~= -1 & (cost < min_cost | min_cost == -1)
           min_cost = cost;
           min_cycle = cycle;
       end
   end

   return

end