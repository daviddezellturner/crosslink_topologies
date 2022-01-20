function datarates = distsToRates(graph)

    datarates = zeros(height(graph),width(graph));

    for i = 1:height(graph)
        for j = 1:width(graph)
            if graph(i,j) == 0
                datarates(i,j) = 0;
            elseif graph(i,j) <= 5
                datarates(i,j) = 50;
            elseif graph(i,j) <= 300
                datarates(i,j) = 10;
            elseif graph(i,j) <= 6000
                datarates(i,j) = 1;
            else
                datarates(i,j) = 0.125;
            end
        end
    end

end

% Slant Range (km)	Data Rate (kbps)
% 0 < x </=5	50
% 5 < x </=300	10
% 300 < x </= 6000	1
% x > 6000	0.125
