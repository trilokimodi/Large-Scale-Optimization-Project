
%%
clc
clearvars
p10
com = sortrows(com, 1);
num_nodes = dimX * dimY * 2;
pi = zeros(dimX*dimY*2, 1);
lagrangian_multiplier = ones(dimX*dimY*2, 1);
lagrangian_multiplier = lagrangian_multiplier * (1 / num_nodes);
stepsizeoffset = 1.99999;
h_history = [];

ergodics = zeros(num_nodes, 1);
ergodicscom = com;
pov = 0;
bestpov = 0;
atiteration = 0;

for iteration = 1:1000
    % stepsize offset change after 10 iterations
    if (mod(iteration, 10) == 0)
        stepsizeoffset = 0.95 * stepsizeoffset;
    end
    
    % initializing the values of pi = lagrangian multiplier for every iteration
    pi = lagrangian_multiplier;
    
    nl = gsp(dimX, dimY, pi, k, com);
    
    last = 0;
    okcom = [];
    newnl = [];
    for i = 1:k
        first = last + 1;
        slask = find(nl(last+1:length(nl)) == com(i, 1));
        last = slask(1) + first - 1;
        if (sum(pi(nl(first:last))) < 1)
            okcom = [okcom, i]; %#ok<AGROW>
            newnl = [newnl; nl(first:last)]; %#ok<AGROW>
        end
    end
    
    % calculation of h
    forh = zeros(k, 1);
    j = 1;
    for i = 1:k
        if (j <= numel(newnl))
            while (j <= numel(newnl) && okcom(i) ~= newnl(j))
                forh(i) = forh(i) + lagrangian_multiplier(newnl(j));
                j = j + 1;
            end
        end
        forh(i) = 1 - forh(i);
        j = j + 1;
    end
    maxforh = max(forh);
    h = sum(lagrangian_multiplier) + maxforh;
    
    % step 3%
    % adding for all occurrences of node i.e outgoing and incoming
    subgrad = zeros(num_nodes, 1);
    for i = 1:num_nodes
        for j = 1:numel(newnl)
            if (i == newnl(j))
                subgrad(i) = subgrad(i) + 1;
            end
        end
    end
    
    % removing the outgoing occurences
    for i = 1:num_nodes
        for j = 1:(numel(okcom))
            if (i == okcom(j))
                subgrad(i) = subgrad(i) - 1;
            end
        end
    end
    
    % calculating subgrad i.e. 1-summation,summation for all i's
    for i = 1:num_nodes
        subgrad(i) = 1 - subgrad(i);
    end
    
    % calculating total subgradient for step length. subgrad squared and added
    totsubgrad = sum(subgrad .* subgrad);
    
    % calculating step length
    step_length = (stepsizeoffset * h) / totsubgrad;
    
    % calculating lagrangian multiplier
    for i = 1:num_nodes
        lagrangian_multiplier(i) = max(0, lagrangian_multiplier(i)-step_length*subgrad(i));
    end
    
    h_history = [h_history, h]; %#ok<AGROW>
    
    % step 4 - Ergodics
    % calculation of ergodics pi values
%     ergodicscom = com;
%     rule = 4;
%     summation1 = 0;
%     summation2 = 0;
%     if (iteration >= 2)
%         for i = 1:num_nodes
%             for s = 0:iteration - 2
%                 summation1 = summation1 + (s + 1)^rule;
%             end
%             for s = 0:iteration - 1
%                 summation2 = summation2 + (s + 1)^rule;
%             end
%             firstterm = (summation1 / summation2) * ergodics(i);
%             secondterm = ((iteration^rule) / summation2) * lagrangian_multiplier(i);
%             ergodics(i) = firstterm + secondterm;
%         end
%         
%         % best lower bound - Heuristic approach
%         % using ergodics pi value, we find route for all pairs
%         ergodicslist = gsp(dimX, dimY, ergodics, k, com);
%         
%         %Counting the num of nodes(outgoing nodes) each pair is taking up
%         count = zeros(dimX, 1);
%         j = 1;
%         blockposition = zeros(dimX, 1);
%         for i = 1:k
%             if (j <= numel(ergodicslist))
%                 blockposition(com(i)) = j;
%                 while (com(i) ~= ergodicslist(j))
%                     j = j + 1;
%                     count(com(i)) = count(com(i)) + 1;
%                 end
%             end
%             j = j + 1;
%         end
%         %Finding the Contact pair which takes maximum nodes in sequence - We should change
%         %this to contact pair with maxmum common pairs
%         for i = 1:numel(find(count)) > 0
%             [maxcount, maxcountposition] = max(count);
%             j = blockposition(maxcountposition);
%             checkcount = 0;
%             position = find(com == maxcountposition);
%             while (com(position) ~= ergodicslist(j))
%                 for check = 1:numel(ergodicslist)
%                     if (ergodicslist(j) == ergodicslist(check))
%                         checkcount = checkcount + 1;
%                     end
%                 end
%                 j = j + 1;
%             end
%             checkcount = checkcount - count(com(position));
%             count(maxcountposition) = -1 * count(maxcountposition);
%             %Removal of contact pair if there are sharing nodes
%             if (checkcount >= 1)
%                 count(maxcountposition) = 0;
%                 j = blockposition(maxcountposition);
%                 while (com(position) ~= ergodicslist(j))
%                     ergodicslist(j) = -1;
%                     j = j + 1;
%                 end
%                 ergodicslist(j) = -1;
%                 position = find(ergodicscom == maxcountposition);
%                 ergodicscom(position, :) = [];
%             end
%         end
%         ergodicslist(ergodicslist == -1) = [];
%         for i = 1:numel(count)
%             count(i) = -1 * count(i);
%         end
%         
%         %Finding the best primal obj value - unncessary
%         if (pov < sum(count))
%             pov = sum(count);
%             solutionvector = ergodicscom;
%             bestlist = ergodicslist;
%             bestpi = ergodics;
%         end
%         
%         %Finding the max number of contact pairs that fit in all iterations
%         if (iteration == 2)
%             noc = numel(ergodicscom);
%             maxconpov = sum(count);
%         end
%         if (noc <= numel(ergodicscom))
%             if (noc == numel(ergodicscom))
%                 if (sum(count) < maxconpov)
%                     maxconpov = sum(count);
%                     maxconcom = ergodicscom;
%                     maxconlist = ergodicslist;
%                     maxconpi = ergodics;
%                     atiteration = iteration;
%                 end
%             else
%                 noc = numel(ergodicscom);
%                 maxconpov = sum(count);
%                 maxconcom = ergodicscom;
%                 maxconlist = ergodicslist;
%                 maxconpi = ergodics;
%                 atiteration = iteration;
%             end
%         end
%     end
end

shift = 25;
visagrid(dimX, dimY, nl, com, pi, shift);

% include a plot of h values in the bottom
plot(h_history)
%%
visagrid(dimX,dimY,maxconlist,maxconcom,maxconpi,shift);
