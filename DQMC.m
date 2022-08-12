clear;
tic;

rng(10);

global gamma ita lambda L M

% constants
gamma1 = 1+sqrt(6)/3;
gamma2 = 1-sqrt(6)/3;
ita1 = sqrt(2*(3-sqrt(6)));
ita1m = -sqrt(2*(3-sqrt(6)));
ita2 = sqrt(2*(3+sqrt(6)));
ita2m = -sqrt(2*(3+sqrt(6)));
gamma = [gamma1 gamma1 gamma2 gamma2];
ita = [ita1 ita1m ita2 ita2m];

%parameters
N = 1;
dt = 0.05; % 1/8
L = 4; % 6
U = 4; % 4
num_w = 200; % 200
num_m = 1000; % 1000
beta = 10;
M = round(beta/dt);
lambda = sqrt(U*dt/2);
re_cal = 5;
mul_int = 10; % 必须能整除M

% observables
n_store = zeros(num_m,1);
s2_store = zeros(num_m,1);
Ek_store = zeros(num_m,1);
Eint_store = zeros(num_m,1);
G_mean = zeros(L);

p1_store = zeros(M*L*num_m,1);
p2_store = zeros(M*L*num_m,1);
p3_store = zeros(M*L*num_m,1);
p4_store = zeros(M*L*num_m,1);

Tij = zeros(L);
for i = 1:L-1
    Tij(i,i+1) = -1;
    Tij(i+1,i) = -1;
end
Tij(1,L) = -1;
Tij(L,1) = -1;

stat_store = [1 -1 2 -2];
field0 = ceil(rand(L,M)*4);
V0 = ita(field0);
B0 = cell(M,1);
temp = expm(-dt*Tij);
for i = 1:M
    B0{i} = temp*diag(exp(1i*lambda*V0(:,i)));
end

field = field0;
V = V0;
B = B0;

% % 验证配分函数是否为实数
% Z = det(inv(G{1}))^(2*N)*prod(exp(-1i*lambda*ita(field)*N).*gamma(field),'all')

% calculate Green function
% G = gen_Green(B);

G = gen_Green2(B,mul_int);
% diff = G{1} - GG{1}

% warm-up
 for n = 1:num_w
    for p = 1:M
        for i = 1:L
            it = field(i,p);
            it_new = it;
            while it_new == it
                it_new = ceil(rand*4);
            end
            
            p1 = 1-G{p}(i,i);
            d_ita = ita(it_new)-ita(it);
            p2 = exp(1i*lambda*(d_ita))-1;
            p3 = exp(-1i*lambda*N*(d_ita));
            p4 = gamma(it_new)/gamma(it);
            pro = (1+p1*p2)^(2*N)*p3*p4;
            
            % heat-bath probability
            pro = pro/(1+pro);
            
            if abs(pro) > 1e5
                warndlg("diverge!")
                return
            end
            
            if real(pro) > rand
                field(i,p) = it_new;
                V(i,p) = ita(it_new);
                B{p} = temp*diag(exp(1i*lambda*V(:,p)));                               
                G{p} = update_Green(G{p}, i, d_ita);
%                 G_test = gen_Green2(B,mul_int);
%                 diff = G_test{p} - G{p}
            end            
        end
        
        % iteration of next Green function
        if p < M
            G{p+1} = B{p}*G{p}/B{p};
%             G_test = gen_Green2(B,mul_int);
%             diff = G_test{p} - G{p}
        else
            G{1} = B{M}*G{M}/B{M};
%             G_test = gen_Green2(B,mul_int);
%             diff = G_test{1} - G{1}

%             G_test = gen_Green(B);
%             diff = G_test{1} - G{1}
        end
    end
    
    % recalculate the Green function
    if mod(n,re_cal) == 0
%         G = gen_Green(B);
        G = gen_Green2(B,mul_int);
    end
 end

count = 1;
% sampling
 for n = 1:num_m
    for p = 1:M
        for i = 1:L
            it = field(i,p);
            it_new = it;
            while it_new == it
                it_new = ceil(rand*4);
            end
            
            p1 = 1-G{p}(i,i);
            d_ita = ita(it_new)-ita(it);
            p2 = exp(1i*lambda*(d_ita))-1;
            p3 = exp(-1i*lambda*N*(d_ita));
            p4 = gamma(it_new)/gamma(it);
            pro = (1+p1*p2)^(2*N)*p3*p4;
            
            p1_store(count) = p1;
            p2_store(count) = p2;
            p3_store(count) = p3;
            p4_store(count) = p4;
            count = count +1;
            
            % heat-bath probability
            pro = pro/(1+pro);
            
            if abs(pro) > 1e3
                warndlg("diverge!")
                return
            end
            
            if real(pro) > rand
                field(i,p) = it_new;
                V(i,p) = ita(it_new);
                B{p} = temp*diag(exp(1i*lambda*V(:,p)));                               
                G{p} = update_Green(G{p}, i, d_ita);
%                 G_test = gen_Green(B);
%                 diff = G_test{p} - G{p}
            end
            
        end
        
        % iteration of next Green function
        if p < M
            G{p+1} = B{p}*G{p}/B{p};
        else
            G{1} = B{M}*G{M}/B{M};
%             G_test = gen_Green(B);
%             diff = G_test{1} - G{1}
        end
    end
    
    % recalculate the Green function
    if mod(n,re_cal) == 0
%         G = gen_Green(B);
        G = gen_Green2(B,mul_int);
    end
    
    G_it = G{1};
    n_store(n) = 1 - real(sum(diag(G{1})))/L;
    G_diag = diag(G_it);
    s2_store(n) = (1-2*mean((1-G_diag).^2))*3/4;
    G_mean = G_mean + G_it;
    
    Ek = 0;
    
    for k = 1:L-1
        Ek = Ek + G_it(k,k+1) + G_it(k+1,k);
    end
    Ek = Ek + G_it(L,1) + G_it(1,L);
    Ek_store(n) = 2*N/L*Ek;
    
    Eint_store(n) = U*N*(2*N-1)*mean(real(G_diag) - abs(G_diag).^2 - 1/2);
 end
 
n_mean = mean(n_store)
Ek_mean = mean(Ek_store)
Eint_mean = mean(Eint_store)
s2_mean = mean(s2_store)
G_mean = G_mean/num_m;

toc;

% calculate Green function
function G = gen_Green(B)
global M L
G = cell(M,1);
for p = 1:M    
    count = 1;
    dd = B{p};
    for k = p+1:M
        dd = B{k}*dd;
        count = count +1;
    end
    for k = 1:p-1
        dd = B{k}*dd;
        count = count +1;
    end
    G{p} = inv(eye(L)+dd);
end
end

% update Green function
function G_new = update_Green(x, i, d_ita)
global lambda L
G_it = x;
G_new = G_it;
fact = exp(1i*lambda*d_ita) - 1;
den = 1+(1-G_it(i,i))*fact;
for j = 1:L
    for k = 1:L
        if i == j
            G_new(j,k) = G_new(j,k) - ((1-G_it(j,i))*fact*G_it(i,k))/den;
        else
            G_new(j,k) = G_new(j,k) + (G_it(j,i)*fact*G_it(i,k))/den;
        end
    end
end
end

function G = gen_Green2(B,mul_int)
global M L
G = cell(M,1);
BB = cell(M,1);
cy = M/mul_int;
for p = 1:M
    for k = p:M
        BB{k-p+1} = B{k};
    end
    for k = 1:p-1
        BB{M-p+1+k} = B{k};
    end
    
    count = 1;
    dd = BB{1};
    for k = 2:mul_int
        dd = BB{k}*dd;
        count = count +1;
    end
    [U0,R0] = qr(dd);
    D0 = diag(diag(R0));
    R0 = D0\R0;
    temp = U0;
    
    for i = 2:cy
        d2 = BB{(i-1)*mul_int+1};
        count = count +1;
        for k = (i-1)*mul_int+2:i*mul_int
            d2 = BB{k}*d2;
            count = count +1;
        end
        [U,R] = qr(d2);
        D = diag(diag(R));
        R = D\R;
        temp = R*temp;
        temp = D*temp;
        temp = U*temp;
    end
    temp = temp*D0;
    temp = temp*R0;
    
%     [Q1, R1] = qr(A1);
%     [Q2, R2] = qr(A2);
%     Q2_inv = inv(Q2);
% %     Q1_inv = inv(Q1);
% %     temp = Q2_inv*Q1_inv + R2*R1;
% %     G{p} = Q1_inv*inv(temp)*Q2_inv;
%     temp = Q2_inv/Q1 + R2*R1;
%     G{p} = Q2\inv(temp)/Q1;
    G{p} = inv(eye(L)+temp);
end
end

