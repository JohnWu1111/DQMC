clear;
format long
tic;

U = 4;
L = 4;
LL = L^2;
beta = 0.1;
len = 2^LL;

Tij = sparse(len,len);
I = sparse(eye(2));
K = sparse([1 0; 0 -1]);
C_dag = sparse([0 0; 1 0]);
C = sparse([0 1; 0 0]);
I_len = diag(sparse(ones(len,1)));
I_len2 = diag(sparse(ones(len^2,1)));

count = 0;
% 横向
% pos=1单独赋值
temp = C_dag; 
temp = kron(temp,C);
count = count + 1;
for j = 3:LL
    temp = kron(temp,I);    
end
Tij = Tij - (temp + temp');

% pos=1,...,L-1单独赋值
for i = 2:L-1
    temp = I;
    for j = 2:i-1
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag);
    temp = kron(temp,C);
    count = count + 1;
    for j = i+2:LL
        temp = kron(temp,I);
    end
    Tij = Tij - (temp + temp');
end

% pos=L单独赋值
temp = C_dag; 
for j = 2:L-1
    temp = kron(temp,K);
end
temp = kron(temp,C);
for j = L+1:LL
    temp = kron(temp,I);
end
count = count + 1;
Tij = Tij - (temp + temp');

for i = 2:L
    for j = 1:L-1
        pos = (i-1)*L+j;
        temp = I;
        for k = 2:pos-1
            temp = kron(temp,I);
        end
        temp = kron(temp,C_dag);
        temp = kron(temp,C);
        count = count + 1;
        for k = pos+2:LL
            temp = kron(temp,I);
        end
        Tij = Tij - (temp + temp');
    end
    
    % 右边缘单独赋值
    pos = i*L;
    temp = I;
    for k = 2:pos-L
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag);
    for k = pos-L+2:pos-1
        temp = kron(temp,K);
    end
    temp = kron(temp,C);
    count = count + 1;
    for k = pos+1:LL
        temp = kron(temp,I);
    end
    Tij = Tij - (temp + temp');
end

% 竖向
% pos=1单独赋值
temp = C_dag;
for j = 2:L
    temp = kron(temp,K);
end
temp = kron(temp,C);
count = count + 1;
for j = L+2:LL
    temp = kron(temp,I);
end
Tij = Tij - (temp + temp');

for i = 2:LL-L
    temp = I;
    for j = 2:i-1
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag);
    for j = i+1:i+L-1
        temp = kron(temp,K);
    end
    temp = kron(temp,C);
    count = count + 1;
    for j = i+L+1:LL
        temp = kron(temp,I);
    end
    Tij = Tij - (temp + temp');
end

% pos=LL-L+1单独赋值
temp = C_dag; 
for j = 2:LL-L
    temp = kron(temp,K);
end
temp = kron(temp,C);
count = count + 1;
for j = LL-L+2:LL
    temp = kron(temp,I);
end
Tij = Tij - (temp + temp');

% pos>LL-L+1
for i = LL-L+2:LL
    temp = I;
    for j = 2:(i-LL+L)-1
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag);
    for j = (i-LL+L)+1:i-1
        temp = kron(temp,K);
    end
    temp = kron(temp,C);
    count = count + 1;
    for j = i+1:LL
        temp = kron(temp,I);
    end
    Tij = Tij - (temp + temp');
end

% TT = -2*(kron(kron(kron(C_dag,C),I),I) + kron(kron(kron(C_dag,K),C),I) + kron(kron(kron(I,C_dag),K),C) + kron(kron(kron(I,I),C_dag),C));
% TT = TT + TT';

H1 = kron(Tij,I_len) + kron(I_len,Tij);

% 格点的粒子数算符
temp_N = cell(LL,1);

% 测试量：n上-n下
test_M = cell(LL,1);

% pos=1单独赋值
H2 = sparse(len^2,len^2);
temp = C_dag*C;
for j = 2:LL
    temp = kron(temp,I);
end
temp1 = kron(temp,I_len) + kron(I_len,temp) - I_len2;
temp1 = temp1^2;
temp_N{1} = kron(temp,I_len) + kron(I_len,temp);
test_M{1} = (kron(temp,I_len) - kron(I_len,temp))^2;
H2 = H2 + temp1*U/2;

for i = 2:LL
    temp = I;
    for j = 2:i-1
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag*C);
    for j = i+1:LL
        temp = kron(temp,I);
    end
    
    temp1 = kron(temp,I_len) + kron(I_len,temp) - I_len2;
    temp1 = temp1^2;
    temp_N{i} = kron(temp,I_len) + kron(I_len,temp);
    test_M{i} = (kron(temp,I_len) - kron(I_len,temp))^2 ;
    H2 = H2 + temp1*U/2;
end

H = H1 + H2;

[V,e] = eigs(H,1,'smallestreal');

% e = diag(D);
% 
% z = sum(exp(-beta*e));

test_n = zeros(1,LL);

for j = 1:LL
    test_it = test_M{j};
    test_n(j) = diag(V'*test_it*V);
end
E_k_n = diag(V'*H1*V)/LL;
E_int_n = diag(V'*H2*V)/LL - U/2;

% E_k = sum(E_k_n.*exp(-beta*e))/z;
% E_int = sum(E_int_n.*exp(-beta*e))/z;
% N = sum(N_n.*exp(-beta*e))/z;
% test = sum(test_n.*exp(-beta*e))/z;
test_f = mean(test_n)*3/4;

G_M = cell(LL);
% pos=1单独赋值 
for j = 2:LL
    temp = C_dag;
    for k = 2:j-1
        temp = kron(temp,K);
    end
    temp = kron(temp,C);
    for k = j+1:LL
        temp = kron(temp,I);
    end
    G_M{1,j} = kron(temp,I_len);
end

for i = 2:LL-1
    for j = i+1:LL
        temp = I;
        for k = 2:i-1
            temp = kron(temp,I);
        end
        temp = kron(temp,C_dag);    
        for k = i+1:j-1
            temp = kron(temp,K);
        end
        temp = kron(temp,C);
        for k = j+1:LL
            temp = kron(temp,I);
        end
        G_M{i,j} = kron(temp,I_len);
    end
end

% pos=1单独赋值 
G_M{1,1} = kron(C_dag*C,diag(sparse(ones(4^LL/2,1))));

for i = 2:LL
    temp = I;
    for k = 2:i-1
        temp = kron(temp,I);
    end
    temp = kron(temp,C_dag*C);
    for k = i+1:LL
        temp = kron(temp,I);
    end
    G_M{i,i} = kron(temp,I_len);
end

G = zeros(LL);
G_n = cell(LL);
for i = 1:LL
    for j = i:LL
        G_M_it = G_M{i,j};
        G_n{i,j} = diag(V'*G_M_it*V);
%         G(i,j) = sum(G_n{i,j}.*exp(-beta*e))/z;
        G(i,j) = G_n{i,j};
        if i ~= j
            G(j,i) = conj(G(i,j));
        end
    end
end

Ek = 0;
for k = 1:LL-1
    Ek = Ek + G(k,k+1) + G(k+1,k);
end
Ek = Ek + G(LL,1) + G(1,LL);

Ek = -2/LL*Ek;

toc;
