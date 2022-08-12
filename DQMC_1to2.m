clear;
tic;
format long
warning('off');

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
dt = 0.02; % 1/8
L = 4; % 6
U = 4; % 4
num_w = 200; % 200
num_m = 1000; % 1000
beta = 10;
M = round(beta/dt);
lambda = sqrt(U*dt/2);
re_cal = 20;
mul_int = 10; % 必须能整除M
output = num_m/100;

% observables
n_store = zeros(num_m,1);
s2_store = zeros(num_m,1);
Ek_store = zeros(num_m,1);
Eint_store = zeros(num_m,1);
G_mean = zeros(L);

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
eil = exp(1i*lambda);
eV0 = eil.^(V0');
B0 = cell(M,1);
temp = expm(-dt*Tij);
tempm = expm(-dt*Tij/2);
tempp = expm(dt*Tij/2);
for p = 1:M
    B0{p} = temp.*eV0(p,:);
end

field = field0;
V = V0;
eV = eV0;
B = B0;

% % 验证配分函数是否为实数
% Z = det(inv(G{1}))^(2*N)*prod(exp(-1i*lambda*ita(field)*N).*gamma(field),'all')

% calculate Green function
G = cell(M,1);
G{1} = gen_Green(B,1);

% G{1} = gen_Green3(B,mul_int,1);
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
            %             p2 = exp(1i*lambda*(d_ita))-1;
            %             p3 = exp(-1i*lambda*N*(d_ita));
            p2 = eil^d_ita-1;
            p3 = 1/eil^(N*d_ita);
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
                eV(p,i) = eil^V(i,p);
                B{p} = temp.*eV(p,:);
                G{p} = update_Green(G{p}, i, p2, L);
                %                 G_test = gen_Green2(B,mul_int);
                %                 diff = G_test{p} - G{p}
            end
        end
        
        % iteration of next Green function
        if p < M
            % recalculate the Green function
            if mod(p,re_cal) == 0
                G{p+1} = gen_Green(B,p+1);
                %         G{1} = gen_Green3(B,mul_int,1);
            else
                G{p+1} = B{p}*G{p}/B{p};
            %             G_test = period_Green(G{p},B{p});
            %             G_test = gen_Green2(B,mul_int);
            %             diff = G_test - G{p+1}
            end
        else
            G{1} = B{M}*G{M}/B{M};
            %             G{1} = period_Green(G{M},B{M});
            %             G_test = gen_Green(B);
            %             diff = G_test{1} - G{1}
        end
    end
    
end

count = 1;
% sampling
for n = 1:num_m
    if mod(n,output) == 0
        n/output
    end
    for p = 1:M
        for i = 1:L
            it = field(i,p);
            it_new = it;
            while it_new == it
                it_new = ceil(rand*4);
            end
            
            p1 = 1-G{p}(i,i);
            d_ita = ita(it_new)-ita(it);
            %             p2 = exp(1i*lambda*(d_ita))-1;
            %             p3 = exp(-1i*lambda*N*(d_ita));
            p2 = eil^d_ita-1;
            p3 = 1/eil^(N*d_ita);
            p4 = gamma(it_new)/gamma(it);
            pro = (1+p1*p2)^(2*N)*p3*p4;
            
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
                eV(p,i) = eil^V(i,p);
                B{p} = temp.*eV(p,:);
                G{p} = update_Green(G{p}, i, p2, L);
                %                 G_test = gen_Green(B);
                %                 diff = G_test{p} - G{p}
            end
            
        end
        
        % iteration of next Green function
        if p < M
            % recalculate the Green function
            if mod(p,re_cal) == 0
                G{p+1} = gen_Green(B,p+1);
                %         G{1} = gen_Green3(B,mul_int,1);
            else
                G{p+1} = B{p}*G{p}/B{p};
            %             G_test = period_Green(G{p},B{p});
            %             G_test = gen_Green2(B,mul_int);
            %             diff = G_test - G{p+1}
            end
        else
            G{1} = B{M}*G{M}/B{M};
            %             G{1} = period_Green(G{M},B{M});
            %             G_test = gen_Green(B);
            %             diff = G_test{1} - G{1}
        end
    end
    
    for p = 1:M
        G_it = tempp*G{p}*tempm;
        G_diag = diag(G_it);
        
        n_store(n) = n_store(n) + 1 - real(sum(diag(G{p})))/L;
        
        Ek = 0;
        for k = 1:L-1
            Ek = Ek + G_it(k,k+1) + G_it(k+1,k);
        end
        Ek = Ek + G_it(L,1) + G_it(1,L);
        Ek_store(n) = Ek_store(n) + 2*N/L*Ek;
        
        Eint_store(n) = Eint_store(n) + U*N*(2*N-1)*mean(real(G_diag) - abs(G_diag).^2 - 1/2);
        
        s2_store(n) = s2_store(n) + (1-2*mean((1-G_diag).^2))*3/4;
        
        G_mean = G_mean + G_it;
    end
    
    n_store(n) = n_store(n)/M;
    Ek_store(n) = Ek_store(n)/M;
    Eint_store(n) = Eint_store(n)/M;
    s2_store(n) = s2_store(n)/M;    
    
end

n_mean = mean(n_store)
Ek_mean = mean(Ek_store)
Eint_mean = mean(Eint_store)
s2_mean = mean(s2_store)
G_mean = G_mean/(num_m*M);

G1 = mean(diag(G_mean))
G2 = (G_mean(1,2) + G_mean(2,3) + G_mean(3,4) + G_mean(4,1) + G_mean(2,1) + G_mean(3,2) + G_mean(4,3) + G_mean(1,4))/8
G3 = (G_mean(1,3) + G_mean(2,4) + G_mean(3,1) + G_mean(4,2))/4

toc;

% calculate Green function
function G = gen_Green(B,p)
global M L
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
G = inv(eye(L)+dd);
end

% update Green function
function G_new = update_Green(G_it, i, p2,L)
G_new = G_it;
% fact = exp(1i*lambda*d_ita) - 1;
den = 1+(1-G_it(i,i))*p2;
% for j = 1:L
%     for k = 1:L
%         if i == j
%             G_new(j,k) = G_new(j,k) - ((1-G_it(j,i))*fact*G_it(i,k))/den;
%         else
%             G_new(j,k) = G_new(j,k) + (G_it(j,i)*fact*G_it(i,k))/den;
%         end
%     end
% end
temp = zeros(L,1);
temp(i) = 1;
G_new = G_new + (G_it(:,i)-temp)*G_it(i,:)*p2/den;
% G_new = G_new + G_it(:,i)*G_it(i,:)*fact/den;
% G_new(i,:) = G_new(i,:) - G_it(i,:)*fact/den;
end

function G = gen_Green3(B,mul_int,p)
global M L
cy = M/mul_int;
count = 1;
dd = B{p};
for k = 2:mul_int
    pos = mod(p+k-2,M)+1;
    dd = B{pos}*dd;
    count = count +1;
end
[U0,R0] = qr(dd);
D0 = diag(diag(R0));
R0 = D0\R0;
temp = U0;

for i = 2:cy
    pos = mod(p+(i-1)*mul_int-1,M)+1;
    d2 = B{pos};
    count = count +1;
    for k = (i-1)*mul_int+2:i*mul_int
        pos = mod(p+k-2,M)+1;
        d2 = B{pos}*d2;
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
G = inv(eye(L)+temp);
end

function G_new = period_Green(x,y)
[U1,R1] = qr(x);
D1 = diag(diag(R1));
R1 = D1\R1;
[U2,R2] = qr(y');

R2 = R2';
U2 = U2';
temp1 = U2*U1;
temp2 = R1/U2;

% UDR
D2 = diag(diag(R2'));
R2_o = R2;
R2 = R2/D2;
% diff = y - R2*D2*U2
G_new = R2_o*temp1*D1*temp2/(D2*R2);

% UR
% G_new = R2*temp1*D1*temp2/R2;
end