function [V,N] = UniformVector(varargin)
    N = varargin{1};
    M = varargin{2};
    seed = varargin{3};
    WVMethod = varargin{4};
    if nargin > 4
        data_set = varargin{5};
    end
    switch WVMethod
        case 'UNV'
            mu = zeros(1,M);
            sigma = eye(M,M);
            rng(seed);
            R = mvnrnd(mu,sigma,N);
            V = abs(R./repmat(sqrt(sum(R.^2,2)),1,M));
        case 'DAS'
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'JAS'
            V = zeros(N,M);
            rng(seed);
            V(:,1) = rand(N,1);
            V(:,1) = 1-(V(:,1).^(1/(M-1)));
            for obj = 2:M-1
                V(:,obj) = (1-sum(V(:,1:obj-1),2)).*(1-rand(N,1).^(1/(M-obj)));
            end
            V(:,M) = 1-sum(V,2);
            V = V./sqrt(sum(V.^2,2));
        case 'MSS-U'
            mu = zeros(1,M);
            sigma = eye(M,M);
            rng(seed);
            if M==3
                R = mvnrnd(mu,sigma,49770);
            elseif M==5
                R = mvnrnd(mu,sigma,46376); 
            elseif M==10
                R = mvnrnd(mu,sigma,48620); 
            end
            S = abs(R./sqrt(sum(R.^2,2)));
            V = zeros(N,M);
            V(1:M,:) = eye(M,M);
            for i = M+1:N
               distance = dist(V(1:i-1,:),S');
               distance = min(distance);
               [~,maxInd] = max(distance);
               V(i,:) = S(maxInd,:);
               S(maxInd,:) = [];
            end
        case 'MSS-D'
            H1 = 1;
            ZN = 0;
            if M==3
                ZN = 49770;
            elseif M==5
                ZN = 46376; 
            elseif M==10
                ZN = 48620; 
            end
            while nchoosek(H1+M,M-1) <= ZN
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            if H1 < M
                H2 = 0;
                while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
                    H2 = H2 + 1;
                end
                if H2 > 0
                    W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
                    W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
                    V  = [V;W2/2+1/(2*M)];
                end
            end
            S = V./sqrt(sum(V.^2,2));
            V = zeros(N,M);
            V(1:M,:) = eye(M,M);
            for i = M+1:N
               distance = dist(V(1:i-1,:),S');
               distance = min(distance);
               [~,maxInd] = max(distance);
               V(i,:) = S(maxInd,:);
               S(maxInd,:) = [];
            end
        case 'UNVS' %Symmetric UNV
            mu = zeros(1,M);
            sigma = eye(M,M);
            rng(seed);
            R = mvnrnd(mu,sigma,floor(N/M));
            Vt = abs(R./repmat(sqrt(sum(R.^2,2)),1,M));
            V = zeros(floor(N/M)*M,M);
            for i = 1:M
                V((i-1)*floor(N/M)+1:i*floor(N/M),:) = Vt;
                Vt = circshift(Vt,[0,1]);
            end
        case 'JASS' %Symmetric JAS
            Vt = UniformVector(floor(N/M),M,seed,'JAS');
            V = zeros(floor(N/M)*M,M);
            for i = 1:M
                V((i-1)*floor(N/M)+1:i*floor(N/M),:) = Vt;
                Vt = circshift(Vt,[0,1]);
            end
        case 'S_energy' %S_energy s=small
            data = load(['S-energy_seed-',num2str(seed),'.mat']);
            DV = data.DV;
            if M == 3
                V = DV{1,1}(:,:,1);
            elseif M == 5
                V = DV{1,2}(:,:,1);
            elseif M == 10
                V = DV{1,3}(:,:,1);
            elseif M == 15
                V = DV{1,4}(:,:,1);
            end
            N = size(V,1);
        case 'UDH' %UDH
            prime_max = log(M-2)/(M-2);
            num_prime = length(primes(prime_max));
            while(num_prime ~= M-2)
                if num_prime < M-2
                    prime_max = prime_max+1;
                elseif num_prime > M-2
                    prime_max = prime_max-1;
                end
                num_prime = length(primes(prime_max));
            end
            p = primes(prime_max);
            V = zeros(N,M);
            for i = 1:N
                V(i,1) = (2*i-1)/(2*N);
                for j = 2:M-1
                    f = 1/p(1,j-1);
                    d = i;
                    while d>0
                        V(i,j) = V(i,j)+f*(mod(d,p(1,j-1)));
                        d = floor(d/p(1,j-1));
                        f = f/p(1,j-1);
                    end
                end
            end
            Vt = V;
            for t = 1:N
                for i = 1:M-1
                    dot = 1;
                    for j = 1:i-1
                        dot = dot*V(t,j)^(1/(M-j));
                    end
                    Vt(t,i) = (1-V(t,i)^(1/(M-i)))*dot;
                end
                Vt(t,M) = 1;
                for j = 1:M-1
                    Vt(t,M) = dot*V(t,j)^(1/(M-j));
                end
            end
            V = Vt;
            V = abs(V./repmat(sqrt(sum(V.^2,2)),1,M));
        case 'CDAS' %True Center DAS method
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            H1 = H1+1;
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V(sum(V==H1,2)>0,:) = [];
            V = V+0.5;
            V = ([V,zeros(size(V,1),1)+H1-0.5]-[zeros(size(V,1),1),V(:,1:M-1)-0.5])/(H1+(M-2)*0.5);
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'DASNZ'
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            V = max(V,1e-6);
            N = size(V,1);
        case 'MSS-UNZ'
            mu = zeros(1,M);
            sigma = eye(M,M);
            rng(seed);
            if M==3
                R = mvnrnd(mu,sigma,49770);
            elseif M==5
                R = mvnrnd(mu,sigma,46376); 
            elseif M==10
                R = mvnrnd(mu,sigma,48620); 
            end
            S = abs(R./sqrt(sum(R.^2,2)));
            V = zeros(N,M);
            V(1:M,:) = eye(M,M);
            for i = M+1:N
               distance = dist(V(1:i-1,:),S');
               distance = min(distance);
               [~,maxInd] = max(distance);
               V(i,:) = S(maxInd,:);
               S(maxInd,:) = [];
            end
            V = max(V,1e-6);
        case 'MSS-DNZ'
            H1 = 1;
            ZN = 0;
            if M==3
                ZN = 49770;
            elseif M==5
                ZN = 46376; 
            elseif M==10
                ZN = 48620; 
            end
            while nchoosek(H1+M,M-1) <= ZN
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            if H1 < M
                H2 = 0;
                while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
                    H2 = H2 + 1;
                end
                if H2 > 0
                    W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
                    W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
                    V  = [V;W2/2+1/(2*M)];
                end
            end
            S = V./sqrt(sum(V.^2,2));
            V = zeros(N,M);
            V(1:M,:) = eye(M,M);
            for i = M+1:N
               distance = dist(V(1:i-1,:),S');
               distance = min(distance);
               [~,maxInd] = max(distance);
               V(i,:) = S(maxInd,:);
               S(maxInd,:) = [];
            end
            V = max(V,1e-6);
        case 'Boundary' %boundary of DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            H1 = 40;
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            V = V(sum(V~=0,2)<M,:);
            N = size(V,1);
        case 'Inner' %Inner of DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            H1 = 17;
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            V = V(sum(V~=0,2)==M,:);
            N = size(V,1);
        case 'UCD' %DAS and CDAS
            das  = UniformVector(N,M,seed,'DASPlus');
            cdas = UniformVector(N,M,seed,'CDASPlus');
            das  = das(sum(das~=0,2)<M,:);
            num = size(cdas,1)+size(das,1);
            if M==10
                ratio = 0.8;
            elseif M==5
                ratio = 0.4;
            elseif M==3
                ratio = 0.2;
            end
            num_das = floor(N*ratio);
            num_cdas = N-num_das;
            ind_das = select(das,num_das,seed);
            ind_cdas = select(cdas,num_cdas,seed);
            V = [das(ind_das,:);cdas(ind_cdas,:)];
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'Proposed'%EquArea rand distribution
            if M == 3
                eud = sqrt(1/3+2*(1/sqrt(2)-1/sqrt(3))^2);
                l_max = 2*asin(eud/2);
                d = l_max:-0.00001:0;
                d_sel = zeros(1,N+1);
                count = 1;
                for i = 1:length(d)
                    a = Area(d(1,i));
                    if pi/6-a >= pi*count/(6*N)
                        d_sel(1,count+1) = d(1,i);
                        count = count+1;
                    end
                end
                d_sel(1,1) = l_max;
                V = zeros(N,3);
                rng(seed);
                for i = 1:N
                    d_max = d_sel(1,i);
                    d_min = d_sel(1,i+1);
                    A_max = sqrt(3)-2*sqrt(3)*(sin(d_max/2))^2;
                    A_min = sqrt(3)-2*sqrt(3)*(sin(d_min/2))^2;
                    C_min = (A_max-sqrt(6-2*A_max^2))/3;
                    C_max = (A_min-sqrt(6-2*A_min^2))/3;
                    OT = zeros(M);
                    for j = 1:M
                        count = 1;
                        OT(j,j) = rand*(C_max-C_min)+C_min;
                        C_e = OT(j,j);
                        S_e = sqrt(1-C_e^2);
                        x_max = S_e*sqrt(1-C_e^2/S_e^2);
                        x_min = S_e*C_e/S_e;
                        count = count+1;
                        for k = 1:M
                            if k ~= j && count < M
                                OT(j,k) = rand*(x_max-x_min)+C_min;
                                count = count+1;
                            elseif k ~= j && count == M
                                OT(j,k) = sqrt(1-sum(OT(j,OT(j,:)>0).^2,2));
                            end
                        end
                    end
                    r = floor(rand*3+1);
                    V(i,:) = OT(r,:);
                end
            end
        case 'DASPlus'
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            N = size(V,1);
        case 'CDASPlus' %True Center DAS method
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            H1 = H1+2;
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V(sum(V==H1,2)>0,:) = [];
            V = V+0.5;
            V = ([V,zeros(size(V,1),1)+H1-0.5]-[zeros(size(V,1),1),V(:,1:M-1)-0.5])/(H1+(M-2)*0.5);
            N = size(V,1);
        case 'Proposed2'%EquDist uniform distribution
            if M == 3
                eud = sqrt(1/3+2*(1/sqrt(2)-1/sqrt(3))^2);
                l_max = 2*asin(eud/2);
                d_sel = linspace(l_max,0,N+1);
                V = zeros(N,3);
                rng(seed);
                for i = 1:N
                    d_max = d_sel(1,i);
                    d_min = d_sel(1,i+1);
                    A_max = sqrt(3)-2*sqrt(3)*(sin(d_max/2))^2;
                    A_min = sqrt(3)-2*sqrt(3)*(sin(d_min/2))^2;
                    C_min = (A_max-sqrt(6-2*A_max^2))/3;
                    C_max = (A_min-sqrt(6-2*A_min^2))/3;
                    OT = zeros(M);
                    for j = 1:M
                        count = 1;
                        OT(j,j) = rand*(C_max-C_min)+C_min;
                        C_e = OT(j,j);
                        S_e = sqrt(1-C_e^2);
                        x_max = S_e*sqrt(1-C_e^2/S_e^2);
                        x_min = S_e*C_e/S_e;
                        count = count+1;
                        for k = 1:M
                            if k ~= j && count < M
                                OT(j,k) = rand*(x_max-x_min)+x_min;
                                count = count+1;
                            elseif k ~= j && count == M
                                OT(j,k) = sqrt(1-sum(OT(j,OT(j,:)>0).^2,2));
                            end
                        end
                    end
                    r = floor(rand*3+1);
                    V(i,:) = OT(r,:);
                end
            end    
        case 'Proposed3'%Sphere DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            t1 = linspace(pi/2,0,H1+1);
            A = [];
            for i = 1:H1+1
                angle = zeros(H1+2-i,2);
                angle(:,2) = t1(1,i);
                if i < H1+1
                    t2 = linspace(pi/2,0,H1+2-i);
                    for j = 1:H1+2-i
                        angle(j,1) = t2(1,j);
                    end
                else
                    angle(H1+2-i,1) = 0;
                end
                A = [A;angle];
            end
            V = zeros(size(A,1),3);
            V(:,3) = cos(A(:,2));
            V(:,2) = sin(A(:,2)).*cos(A(:,1));
            V(:,1) = sin(A(:,2)).*sin(A(:,1));
            V = max(V,1e-6);
            N = size(V,1);
        case 'Proposed4'%Centers of Sphere DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            theta = acos(1/sqrt(3));
            t1 = linspace(pi/2-1/H1*pi/2+theta/H1,theta/H1,H1);
            A = [];
            for i = 1:H1
                angle = zeros(H1-i+1,2);
                angle(:,2) = t1(1,i);
                if i < H1
                    K = H1-i+1;
                    t2 = linspace(pi/2-1/K*pi/4,1/K*pi/4,K);
                    for j = 1:H1-i+1
                        angle(j,1) = t2(1,j);
                    end
                else
                    angle(H1+1-i,1) = pi/4;
                end
                A = [A;angle];
            end
            V = zeros(size(A,1),3);
            V(:,3) = cos(A(:,2));
            V(:,2) = sin(A(:,2)).*cos(A(:,1));
            V(:,1) = sin(A(:,2)).*sin(A(:,1));
            V = max(V,1e-6);
            N = size(V,1);
        case 'Proposed5' %Center of DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            t1 = linspace(1/(3*H1),1-2/(3*H1),H1);
            A = [];
            for i = 1:H1
                angle = zeros(H1-i+1,3);
                angle(:,3) = t1(1,i);
                if i < H1
                    t2 = linspace(1/(3*H1),1-(i-1)/H1-2/(3*H1),H1-i+1);
                    for j = 1:H1-i+1
                        angle(j,1) = t2(1,j);
                        angle(j,2) = 1-angle(j,1)-angle(j,3);
                    end
                else
                    angle(H1+1-i,1) = 1/(3*H1);
                    angle(H1+1-i,2) = 1/(3*H1);
                end
                A = [A;angle];
            end
            V = A;
            V = max(V,1e-6);
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'Proposed6' %Dense Center of DAS
            H1 = 1;
            while nchoosek(H1+M,M-1) <= N
                H1 = H1 + 1;
            end
            t1j = linspace(1/(3*H1),1-2/(3*H1),H1);
            t1o = linspace(2/(3*H1),1-1/(3*H1),H1);
            A = [];
            for i = 1:H1
                angle = zeros(2*(H1-i)+1,3);
                if i < H1
                    t2j = linspace(1/(3*H1),1-(i-1)/H1-2/(3*H1),H1-i+1);
                    t2o = linspace(2/(3*H1),1-1/(3*H1)-i/H1,H1-i);
                    for j = 1:2*(H1-i)+1
                        if mod(j,2) == 1
                            angle(j,3) = t1j(1,i);
                            angle(j,1) = t2j(1,ceil(j/2));
                        else
                            angle(j,3) = t1o(1,i);
                            angle(j,1) = t2o(1,ceil(j/2));
                        end
                        angle(j,2) = 1-angle(j,1)-angle(j,3);
                    end
                else
                    angle(H1+1-i,1) = 1/(3*H1);
                    angle(H1+1-i,2) = 1/(3*H1);
                    angle(H1+1-i,3) = 1-angle(H1+1-i,1)-angle(H1+1-i,2);
                end
                A = [A;angle];
            end
            V = A;
            V = max(V,1e-6);
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'Proposed7' %IGD based direction vector generation method
            H1 = 1;
            ZN = 0;
            if M==3
                ZN = 990;
            elseif M==5
                ZN = 1000; 
            end
            while nchoosek(H1+M,M-1) <= ZN
                H1 = H1 + 1;
            end
            V = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
            V = ([V,zeros(size(V,1),1)+H1]-[zeros(size(V,1),1),V])/H1;
            if H1 < M
                H2 = 0;
                while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
                    H2 = H2 + 1;
                end
                if H2 > 0
                    W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
                    W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
                    V  = [V;W2/2+1/(2*M)];
                end
            end
            V = max(V,1e-6);
            S = V./sqrt(sum(V.^2,2));
            V = S;
            for i = 1:N
                Score = zeros(1,size(S,1));
                for j = 1:size(S,1)
                    VT = [V;S(j,:)];
                    ST = [S(1:j-1,:);S(j+1:size(S,1),:)];
                    Score(1,j) = IGD(VT,ST);
                end
                min_value = min(Score);
                min_index = find(Score==min_value);
                V = [V;S(min_index(1),:)];
                S(min_index(1),:) = [];
            end
        case 'Proposed8' %The gradient of the length of direction vectors
            mu = zeros(1,M);
            sigma = eye(M,M);
            rng(seed);
            if M==3
                R = mvnrnd(mu,sigma,49770);
            elseif M==5
                R = mvnrnd(mu,sigma,46376); 
            end
            S = abs(R./sqrt(sum(R.^2,2)));
            V = zeros(N,M);
            Step = pi/2/N;                 
        case 'Proposed9' %data_set dependent DV method
            V = data_set;
            V = V./sqrt(sum(V.^2,2));
            V = V(1:N,:);
        case 'Proposed10' %True IGD based
            if M == 3
                V = load(['SMSEMOA_DTLZ2_M3_D12_',num2str(seed),'.mat']);
            elseif M == 5
                V = load(['SMSEMOA_DTLZ2_M5_D14_',num2str(seed),'.mat']);
            end
            V = V.result{1,2}.objs;
            V = V./repmat(sqrt(sum(V.^2,2)),1,M);
            N = size(V,1);
        case 'Proposed14' %S_energy s=large
            V = readmatrix(['M-',num2str(M),'_seed-',num2str(seed),'_N-',num2str(N),'_s-',num2str((M-1)*10),'.dat']);
            V = V(:,1:M);
    end
end

function o = norm(val_min,val_max)
    o = normrnd((val_min+val_max)/2,(val_max-val_min)/2);
    while o >= val_max || o <= val_min
        o = normrnd((val_min+val_max)/2,(val_max-val_min)/2);
    end
end
function ind = select(dvs,num,seed)
    H = 1;
    [N,M] = size(dvs);
    while nchoosek(H+M,M-1) <= N
        H = H + 1;
    end
    H = H+1;
    ref = (H+1)/H;
    V = dvs;
    HV = prod(ref-V,2);
    HV = HV/sum(HV);
    ruler = cumsum(HV);
    rng(seed);
    ind = zeros(1,num);
    is_select = false(1,N+1);
    for i = 1:num
        is_select(end) = false;
        while ~is_select(end)
            u = rand;
            front = 0;
            for j = 1:length(ruler)
                if u >= front && u <= ruler(j)
                    tind = j;
                    if ~is_select(tind)
                        ind(1,i) = tind;
                        is_select(tind) = true;
                        is_select(end) = true;
                    end
                    break;
                else
                    front = ruler(j);
                end
            end
        end
    end
end
