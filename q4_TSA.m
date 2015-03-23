clear

N = 30;         %Population size
n = 5;          %Number of cities
Pr = 0.65;      %Probability of reproduction
Pm = 0;         %Probability of mutation
simCount = 500; %Number of simulations

% Uncomment below to generate randome journey times
% Random integer between 0 and 50 to represent journey time (cost)
% cost = randi([0 50], n);

cost = [41    36    26    35    12
        27    17    25    44     1
        21    14    34    50    38
        13    10    37    36    44
        4     1     32    16    35]; 

%Remove cost at diagonal   
for i = 1:length(cost)
    cost(i,i)= NaN;
end
clear i

rs = [];    %Contain routes generated
rCs = [];   %Contain cost for each location pair
rC = [];    %Total route cost
v = perms([1:n]); %Generate all possible routes using permutation

for i = 1:N                             %For each route
    route = v(randi([1 length(v)]),:);  %Select a route randomly from v
    rs(i,:) = route;                    %Store selected route
    c = 0;                              %Initiase cost
    for j = 1:n-1
        tmp = cost(route(j),route(j+1));
        c = c + tmp;                    %Increment route cost
        rCs(i,j) = tmp;                 %Store route cost 
    end   
    rC(i) = c;                          %Store total route cost
end
clear tmp c route i v j 

for sim = 1:simCount
    
    %Implement knowledge based multiple inversion
    for i = 1:N %For each route
        route = rs(i,:); %Select route
        routeCost = rCs(i,:); %Select route cost
        routeCostsorted = sort(routeCost,'descend'); %Sort route cost descending order

        %Generate pairs column in sorted cost order
        for j = 1:n-1 
            index = find(routeCost==routeCostsorted(j));
            index = index(1);
            pairs(j,:) = route(index:index+1); 
            routeCost(index)=NaN;
        end
        
        clear index j

        s = 1; %Pair indexing
        check = [];
        index=[];
        route1 = route;

        for k = 1:(n-1)/2
            pair1 = pairs(s,:);
            pair2 = pairs(s+1,:);
            index(1) = find(route==pair1(1));
            index(2) = find(route==pair2(1));
            if index(2) < index(1) %sort ordering
                index([1 2]) = index([2 1]);
            end

            u = ismember([pair1 pair2], check); %Check if existing pairs overlap check
            overlap = any(u(:)== 1);

            if overlap == 0 %If no overlap
                slice = route(index(1)+1:index(2)); %slice out array
                slice([1 end]) = slice([end 1]); %switch elements
                route1(index(1)+1:index(2))= slice;
                check = [check,slice]; %add the slice in check
            end
            s=s+2; %increment to next pairs
        end
        new_rs(i,:)= route1;
    end

    clear i check k index pair1 pair2 pairs route route1 routeCost routeCostsorted slice k u s


    %Update route total cost and idividual costs in each route
    for i = 1:N
        route = new_rs(i,:); %for N times
        c = 0;
        for j = 1:n-1
            tmp = cost(route(j),route(j+1));
            c = c + tmp;
            %new_rCs(i,j) = tmp;
        end   
        new_rC(i) = c;
    end
    clear j c i tmp
    
    %Quick check - knowledge based inversion should improve journey time
    %sum(new_rC) < sum(rC); 
    
    %Natural selection - Roulette wheel method
    for m=1:N
        f(m) = 1./new_rC(m); %The larger the element, the less probable to be selected
    end
    ff = sum(f);
    for m=1:N
        p(m) = f(m)./ff; %Probability of being selected
    end
    cp(1)= p(1);
    for m=2:N
        cp(m) = p(m)+cp(m-1); %cumlative
    end
    %Equivalent to spinning roulette wheel N times.
    R = rand(1,N); %generate barriar of entry randomly
    for m=1:N
        D = cp - R(m); 
        I = find(D > 0); %return index position of the element which exceeds the R threshold
        el(m,:) = new_rs(I(1),:); %Select the first element out of elite group into the next generation
    end

    clear f p R I D m ff cp
    
    %Crossover/breed, el are the parents
    PA = randi(N,1,N+1); 
    childs=[];
    m = 2;
    while m-1 < N+1 %exit while loop only when N child are generated
        r = rand; %probability of breeding
        if (r < Pr) %cross over only when r is smaller than reproduction probability
                    %randomly choose how much to cross over, and which number to stay
            %J = randi(n, 1,randi([1 n-1],1)); 
            a = el(PA(m),:); %select parent 1
            b = el(PA(m-1),:); %select parent 2

            index = randi(n,1); %choose random place to switch
            if n-index==0
                length=0;
            else
                length = randi(n-index, 1);
            end

            b1 = [b(index+length+1:end),b(1:index+length)];
            eltr = a(index:index+length);
            C = b1(~ismember(b1,eltr));

            bk = n-length-index;
            ft = index-1;

            if ft <= 0
                child = [eltr, C];
            elseif bk <= 0
                child = [C, eltr];
            else
                child = [C(bk+1:end), eltr, C(1:bk)];
            end
            childs(m-1,:)=child;
            m=m+1;
        end
    end
    
    %Let's make some mutants!
    for mu = 1:N
        if rand < Pm %if probability is below probability of mutation
            mdex = randperm(n);
            mdex = mdex(1:2);
            pre_mutant = childs(mu,:);
            pre_mutant([mdex(1) mdex(2)]) = pre_mutant([mdex(2) mdex(1)]);
            childs(mu,:) = pre_mutant;
        end
    end
    
    %Update cost for childs
    for i = 1:N
        route = childs(i,:); %for N times
        c = 0;
        for j = 1:n-1
            tmp = cost(route(j),route(j+1));
            c = c + tmp;
            new_rCs(i,j) = tmp; 
        end   
        new_childC(i) = c; %Update childs cost
    end

rs = childs;
rCs = new_rCs;
totalGenCost(sim) = sum(new_childC); %store generation total cost
averageGenCost(sim) = sum(new_childC)/N; %store avergae generation cost
end

cost
plot(1:sim, averageGenCost);
