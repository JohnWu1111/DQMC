clear;
tic;
format long
warning('off');

rng(11);

global L gamma ita lambda len M

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
len = L^2;
U = 4; % 4
num_w = 100;
num = 200; % 1000
beta = 40;
M = round(beta/dt);
lambda = sqrt(U*dt/2);
re_cal = 10;
mul_int = 10; % 必须能整除M
output = num_w/100;

% observables
n_store = zeros(num,1);
s2_store = zeros(num,1);
Ek_store = zeros(num,1);
Eint_store = zeros(num,1);
G_mean = zeros(len);

Tij = gen_H(1);

stat_store = [1 -1 2 -2];
field0 = ceil(rand(len,M)*4);
V0 = ita(field0);
eil = exp(1i*lambda);
eV0 = eil.^(V0');
B0 = cell(M,1);
% temp = expm(-dt*Tij);
% tempm = expm(-dt*Tij/2);
% tempp = expm(dt*Tij/2);
[V1,D1] = eig(Tij);
temp = V1*diag(exp(-dt*diag(D1)))*V1';
tempm = V1*diag(exp(-dt*diag(D1)/2))*V1';
tempp = V1*diag(exp(dt*diag(D1)/2))*V1';
for p = 1:M
    B0{p} = temp.*eV0(p,:);
end

% for i = 1:M
%     B0{i} = temp*diag(exp(1i*lambda*V0(:,i)));
% end

field = field0;
V = V0;
eV = eV0;
B = B0;

% generate the projection matrix
flux = 1e-3;
s = exp(1i*flux/len);
H = gen_H(s);
[VV,D] = eig(H);

e = real(diag(D));
[e1,I] = sort(e);
VV1 = VV(:,I);

% load('from_fortran.mat')
% VV = VV_f;

PR = VV(:,1:len/2);
% PL = PR';
PR = tempp*PR;
PL = PR';
% PL = PL*tempp;

% calculate Green function
G = cell(M,1);
% [ML, MR] = gen_LR(B,1,PL,PR);
[UL, UR] = gen_LR(B,1,PL,PR,mul_int);
G{1} = eye(len) - UR/(UL*UR)*UL;
% G{1} = gen_Green2(ML,MR);
G{1}(1,1)

% % 验证配分函数是否为实数
% Z = det(ML*MR)^(2*N)*prod(exp(-1i*lambda*ita(field)*N).*gamma(field),'all')

for n = 1:num_w
    if mod(n,output) == 0
        n/output
    end
    for p = 1:M
        for i = 1:len
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
            
%             if abs(pro) > 1e3
%                 warndlg("diverge!")
%                 return
%             end
            pro = abs(pro);
            % heat-bath probability
            pro = pro/(1+pro);            
            
            if real(pro) > rand
                field(i,p) = it_new;
                V(i,p) = ita(it_new);
                eV(p,i) = eil^V(i,p);
%                 B{p} = temp.*exp(1i*lambda*V(:,p)).';
                B{p} = temp.*eV(p,:);
                G{p} = update_Green(G{p}, i, p2, len);
                %                 G_test = gen_Green(B);
                %                 diff = G_test{p} - G{p}
            end
            
        end
        
        % iteration of next Green function
        if p < M
            % recalculate the Green function
            if mod(p,re_cal) == 0
                % G{p+1} = gen_Green(B,p+1);
%                 [ML, MR] = gen_LR(B,p+1,PL,PR);
                [UL, UR] = gen_LR(B,p+1,PL,PR,mul_int);
                G{p+1} = eye(len) - UR/(UL*UR)*UL;
            else
                G{p+1} = B{p}*G{p}/B{p};
            %             G_test = period_Green(G{p},B{p});
            %             G_test = gen_Green2(B,mul_int);
            %             diff = G_test - G{p+1}
            end
        else
%             [ML, MR] = gen_LR(B,1,PL,PR);
            [UL, UR] = gen_LR(B,1,PL,PR,mul_int);
            G{1} = eye(len) - UR/(UL*UR)*UL;
%             G{1} = B{M}*G{M}/B{M};
            %             G{1} = period_Green(G{M},B{M});
            %             G_test = gen_Green(B);
            %             diff = G_test{1} - G{1}
        end
    end  
    
end

for n = 1:num
    for p = 1:M
        for i = 1:len
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
            
%             if abs(pro) > 1e3
%                 warndlg("diverge!")
%                 return
%             end
            
            pro = abs(pro);
            
            % heat-bath probability
            pro = pro/(1+pro);
                       
            if real(pro) > rand
                field(i,p) = it_new;
                V(i,p) = ita(it_new);
                eV(p,i) = eil^V(i,p);
                B{p} = temp.*eV(p,:);
                G{p} = update_Green(G{p}, i, p2, len);
                %                 G_test = gen_Green(B);
                %                 diff = G_test{p} - G{p}
            end
            
        end
        
        % iteration of next Green function
        if p < M
            % recalculate the Green function
            if mod(p,re_cal) == 0
                % G{p+1} = gen_Green(B,p+1);
                %                 [ML, MR] = gen_LR(B,p+1,PL,PR);
                [UL, UR] = gen_LR(B,p+1,PL,PR,mul_int);
                G{p+1} = eye(len) - UR/(UL*UR)*UL;
            else
                G{p+1} = B{p}*G{p}/B{p};
                %             G_test = period_Green(G{p},B{p});
                %             G_test = gen_Green2(B,mul_int);
                %             diff = G_test - G{p+1}
            end
        else
            %             [ML, MR] = gen_LR(B,1,PL,PR);
            [UL, UR] = gen_LR(B,1,PL,PR,mul_int);
            G{1} = eye(len) - UR/(UL*UR)*UL;
            %             G{1} = B{M}*G{M}/B{M};
            %             G{1} = period_Green(G{M},B{M});
            %             G_test = gen_Green(B);
            %             diff = G_test{1} - G{1}
        end
    end
    
    p_min = round(M/4);
    p_max = round(M*3/4);
    p_num = p_max - p_min + 1;
    
    for p = p_min:p_max
        G_it = tempp*G{p}*tempm;
        %         G_it = G{p};
        G_diag = diag(G_it);
        
        n_store(n) = n_store(n) + 1 - real(sum(diag(G{p})))/len;
        
        %         Ek = 0;
        %         for k = 1:L-1
        %             Ek = Ek + G_it(k,k+1) + G_it(k+1,k);
        %         end
        %         Ek = Ek + G_it(L,1) + G_it(1,L);
        %         Ek_store(n) = Ek_store(n) + 2*N/L*Ek;
        
        Ek = 0;
        count = 0;
        for x = 1:L
            for y = x:L
                if Tij(x,y) ~= 0
                    Ek = Ek + G_it(x,y) + G_it(y,x);
                    count = count + 2;
                end
            end
        end
        Ek_store(n) = Ek_store(n) + 2*N/L*Ek;
        
        Eint_store(n) = Eint_store(n) + U*N*(2*N-1)*mean(real(G_diag) - abs(G_diag).^2 - 1/2);
%         Eint_store(n) = Eint_store(n) + U*mean((1-G_diag).^2);
        
        s2_store(n) = s2_store(n) + (1-2*mean((1-G_diag).^2))*3/4;
        
        G_mean = G_mean + G_it;
    end
    
    n_store(n) = n_store(n)/p_num;
    Ek_store(n) = Ek_store(n)/p_num;
    Eint_store(n) = Eint_store(n)/p_num;
    s2_store(n) = s2_store(n)/p_num;
    
end

n_mean = mean(n_store)
Ek_mean = mean(Ek_store)
Eint_mean = mean(Eint_store)
s2_mean = mean(s2_store)
G_mean = G_mean/(num*p_num);

G1 = mean(diag(G_mean))
% G2 = (G_mean(1,2) + G_mean(1,3) + G_mean(2,4) + G_mean(3,4) + G_mean(2,1) + G_mean(3,1) + G_mean(4,2) + G_mean(4,3))/8
% G3 = (G_mean(1,4) + G_mean(2,3) + G_mean(4,1) + G_mean(3,2))/4

G2 = 0;
count = 0;
for i = 1:L
    for j = i:L
        if Tij(i,j) ~= 0
            G2 = G2 + G_mean(i,j) + G_mean(j,i);
            count = count + 2;
        end
    end
end
G2 = G2/count

toc;

function [UL, UR] = gen_LR(B,p,PL,PR,mul_int)
global M
count = 0;
num_L = M-p+1;
cy_L = floor(num_L/mul_int);
[V,D,U] = L_MGS(PL);

for i = 1:cy_L
    pos = M - mul_int*(i-1);
    d2 = B{pos};
    count = count +1;
    for j = pos-1:-1:pos-mul_int+1
        d2 = d2*B{j};
        count = count +1;
    end
    [VV,D,U] = L_MGS(U.*D*d2);
    V = V*VV;
end
if mod(num_L,mul_int)~=0
    pos = M - mul_int*cy_L;
    d2 = B{pos};
    count = count +1;
    for j = pos-1:-1:p
        d2 = d2*B{j};     
        count = count +1;
    end
    [~,~,U] = L_MGS(U.*D*d2);
end
UL = U;

num_R = p-1;
cy_R = floor(num_R/mul_int);
[U,D,V] = R_MGS(PR);
for i = 1:cy_R
    pos = 1 + mul_int*(i-1);
    d2 = B{pos};
    count = count +1;
    for j = pos+1:pos+mul_int-1
        d2 = B{j}*d2;
        count = count +1;
    end
    [U,D,VV] = R_MGS(d2*U.*D');
    V = VV*V;
end
if mod(num_R,mul_int)~=0
    pos = 1 + mul_int*cy_R;
    d2 = B{pos};
    count = count +1;
    for j = pos+1:p-1
        d2 = B{j}*d2;
        count = count +1;
    end
    [U,~,~] = R_MGS(d2*U.*D');
end
UR = U;
end

function [V,D,U] = L_MGS(A)
global len
[Q,R] = qr(A');
U = Q';
U = U(1:len/2,:);
D = diag(R);
V = R'./D;
V = V(:,1:len/2);
end

function [U,D,V] = R_MGS(A)
global len
[U,R] = qr(A);
U = U(:,1:len/2);
D = diag(R);
V = R./D';
V = V(1:len/2,:);
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

function G_new = update_Green2(G_it, i, p2,L)
G_new = G_it;
den = 1+(1-G_it(i,i))*p2;
pha = p2/den;

for k = 1:L
    fac = pha*G_it(i,k);
    for j = 1:L
        delta = 0;
        if j == i
            delta = 1;
        end
        G_new(j,k) = G_new(j,k) - (delta - G_it(j,i)*fac);
    end
end
end

function Tij = gen_H(s)
global L len
Tij = zeros(len);
count = 0;
for i = 1:L-1
    for j = 1:L-1
        pos = (i-1)*L+j;
        Tij(pos,pos+1) = Tij(pos,pos+1)-s;
        Tij(pos+1,pos) = Tij(pos+1,pos)-conj(s);
        Tij(pos,pos+L) = Tij(pos,pos+L)-1;
        Tij(pos+L,pos) = Tij(pos+L,pos)-1;
        count = count +1;
    end
    pos = i*L;
    Tij(pos,pos-L+1) = Tij(pos,pos-L+1)-s;
    Tij(pos-L+1,pos) = Tij(pos-L+1,pos)-conj(s);
    Tij(pos,pos+L) = Tij(pos,pos+L)-1;
    Tij(pos+L,pos) = Tij(pos+L,pos)-1;
    count = count +1;
end
for j = 1:L-1
    pos = (L-1)*L+j;
    Tij(pos,pos+1) = Tij(pos,pos+1)-s;
    Tij(pos+1,pos) = Tij(pos+1,pos)-conj(s);
    Tij(pos,pos+L-len) = Tij(pos,pos+L-len)-1;
    Tij(pos+L-len,pos) = Tij(pos+L-len,pos)-1;
    count = count +1;
end
Tij(len,len-L+1) = Tij(len,len-L+1)-s;
Tij(len-L+1,len) = Tij(len-L+1,len)-conj(s);
Tij(len,L) = Tij(len,L)-1;
Tij(L,len) = Tij(L,len)-1;
count = count +1;
end