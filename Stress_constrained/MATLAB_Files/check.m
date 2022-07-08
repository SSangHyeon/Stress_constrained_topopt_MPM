function [dcn]=check(node,x,dc,rmin)
    dcn = zeros(length(node),1);
    for i = 1:length(node)
            sum = 0;
            sum2 = 0;
            ex = find(abs(node(:,1) - node(i,1)) < rmin);
            ey = find(abs(node(:,2) - node(i,2)) < rmin);
            E = intersect(ex,ey);
            for j = 1:length(E)
                A = [node(i,1), node(i,2)];
                B = [node(E(j),1), node(E(j),2)];
                weight = max(0,rmin - norm(A-B));
                sum = sum + weight;
                sum2 = sum2 + weight * x(E(j)) * dc(E(j));
            end
            dcn(i) = sum2 / (x(i) * sum);
    end