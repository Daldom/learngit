%   Power flow compute.
%   It is the main document.

%% ��������������
function main
% set up the variables
global G B G0 B0 H N J L a b;
global P Q;
global BusNum BranchNum;
global Um Ua;
global Delta;
global Pg Qg Pd Qd;
global PQnum PVnumc DeltaNum;
global Um Ua DeltaUa;
global DeltaMax;
global Jacobi;


COUNTNUMMAX = 100;
mpc = case14;
BusNum = size(mpc.bus,1);
BranchNum = size(mpc.branch,1);


formY%�γɽڵ㵼�ɾ���
VoltageInitial;%�ڵ��ѹ��ʼ��
for CountNum = 1:COUNTNUMMAX %ѭ������
    SetUnbalance%���㲻ƽ����
    GetUnbalanceMax%�������ƽ��ֵ
%�Ƚϲ�ƽ�����뾫��
    %���㾫��Ҫ��
    
    %�����㾫��Ҫ��
        %�ж�С�����õ�������
            FormJacobi%�γ��ſ˱Ⱦ���
            GetRevised%����������̵õ���ѹ����
            GetNewVoltageValue%�õ���ѹ��ֵ
end
%%test

end
%% �γɽڵ㵼�ɾ���
function formY
global G B;
global G0 B0;
global BusNum ;
global BranchNum;

mpc = case14;

%��ʼ���ڵ㵼�ɾ��� BusNum * BusNum
G = zeros(BusNum);
B = zeros(BusNum);
G0 = zeros(BusNum);
B0 = zeros(BusNum);
%��֧·Ϊ��λ����ѭ������
for i = 1:BranchNum
    m = mpc.branch(i,1);
    n = mpc.branch(i,2);
    K = mpc.branch(i,9); 
    if K==0
        K = 1;
    end
    
    Rmn = mpc.branch(i,3);
    Xmn = mpc.branch(i,4);
    Bmn0= mpc.branch(i,5);
    
    G10 = ( 1 - K ) / ( K * K ) * Rmn / ( Rmn * Rmn + Xmn * Xmn );
    B10 =-( 1 - K ) / ( K * K ) * Xmn / ( Rmn * Rmn + Xmn * Xmn );
    G20 = ( K - 1 ) / K * Rmn / ( Rmn * Rmn + Xmn * Xmn );
    B20 =-( K - 1 ) / K * Rmn / ( Rmn * Rmn + Xmn * Xmn );
    Gmn = 1 / K * Rmn / ( Rmn * Rmn + Xmn * Xmn );
    Bmn =-1 / K * Xmn / ( Rmn * Rmn + Xmn * Xmn );
    B00 = Bmn0 / 2;
    
    %��֧·׷�ӷ��γɽڵ㵼�ɾ���
    G(m,m) = G(m,m) + Gmn + G10;
    G(m,n) = - Gmn;
    G(n,m) = - Gmn;
    G(n,n) = G(n,n) + Gmn + G20;
    
    B(m,m) = B(m,m) + Bmn + B10 + B00;
    B(m,n) = - Bmn;
    B(n,m) = - Bmn;
    B(n,n) = B(n,n) + Bmn + B20 + B00;
    
    %�洢�Ե�֧·����
    G0(m,n) = G10;
    G0(n,m) = G20;
    B0(m,n) = B10 + B00;
    B0(n,m) = B20 + B00;
end
end
%% ���ýڵ��ѹ������
function VoltageInitial
global BusNum;
global Um Ua;
global Pg Qg Pd Qd;
global PQnum PVnum;

mpc = case14;
SB  = mpc.baseMVA;
%�������ʼ��
Pg = zeros(BusNum,1);
Qg = zeros(BusNum,1);
Pd = zeros(BusNum,1);
Qd = zeros(BusNum,1);
Um = zeros(BusNum,1);
Ua = zeros(BusNum,1);
PQnum = 0;
PVnum = -1;%����ʱ�Ƽ�ƽ��ڵ㣬�ʼ�ȥ

%�ڵ��ѹ����λ
for i = 1:BusNum 
    inode = mpc.bus(i,1);
    itype = mpc.bus(i,2);
    if itype == 3 || itype == 2
        Um(inode) = mpc.bus(i,8);
        Ua(inode) = mpc.bus(i,9);
        PVnum = PVnum + 1;
    else
        Um(inode) = mpc.bus(i,8);
        Ua(inode) = mpc.bus(i,9);
        PQnum = PQnum + 1;
    end
    %���ɹ���
    Pd(inode) = mpc.bus(i,3) / SB;%��Ϊ����ֵ
    Qd(inode) = mpc.bus(i,4) / SB;
end
%�����
Num = size(mpc.gen,1);
for i = 1:Num
    inode = mpc.gen(i,1);
    Pg(inode) = mpc.gen(i,2) / SB;%��Ϊ����ֵ
    Qg(inode) = mpc.gen(i,3) / SB;
end

end
%% ���㹦�ʲ�ƽ����
function SetUnbalance
global BusNum;
global Um Ua DeltaUa;
global G B;
global Delta;
global P Q Pg Qg Pd Qd;
global PQnum PVnum DeltaNum;

mpc = case14;
P = zeros(BusNum,1);
Q = zeros(BusNum,1);
DeltaP = zeros(BusNum,1);
DeltaQ = zeros(BusNum,1);
DeltaUa = zeros(BusNum);
DeltaNum = PQnum * 2 +  PVnum;
Delta = zeros(DeltaNum,1);
line = 1;%��ʼ����line��ʾ�洢��Delat���������

for i = 2:BusNum
    itype = mpc.bus(i,2);
    %���㹦�ʷ���
    for j = 1:BusNum
         DeltaUa(i,j) = Ua(i) - Ua(j);
         P(i) = P(i) + Um(i) * Um(j) * ( G(i,j) * cosd(DeltaUa(i,j)) + B(i,j) * sind(DeltaUa(i,j)) );
         Q(i) = Q(i) + Um(i) * Um(j) * ( G(i,j) * sind(DeltaUa(i,j)) - B(i,j) * cosd(DeltaUa(i,j)) );
    end
    
    if itype == 1 %PQ�ڵ���������
        DeltaP(i) = Pg(i) - Pd(i) - P(i);
        DeltaQ(i) = Qg(i) - Qd(i) - Q(i); 
        %Delta(line) = DeltaP(i);
        %Delta(line+1) = DeltaQ(i);
        %line = line + 2;
    end    
    if itype == 2
        DeltaP(i) = Pg(i) - Pd(i) - P(i);
        %Delta(line) = DeltaP(i);
        %line = line + 1;
    end
    
end
for i = 2:BusNum
    itype = mpc.bus(i,2);
    if itype == 1 %PQ�ڵ���������
        Delta(line)   = DeltaP(i);
        Delta(line+1) = DeltaQ(i);
        line = line + 2;   
    elseif itype == 2
        Delta(line) = DeltaP(i);
        line = line + 1;
    end
end  
end
%% �������ƽ��ֵ
function GetUnbalanceMax
global DeltaMax;
global Delta;
global DeltaNum;

DeltaMax = 0;
AbsDelta = zeros(DeltaNum,1);
for i = 1:DeltaNum
    if Delta(i)<0
        AbsDelta(i) = -Delta(i);
    else
        AbsDelta(i) =  Delta(i);
    end
    
    if AbsDelta(i)>DeltaMax %�����ƽ����
        DeltaMax = AbsDelta(i);
    end
end
end
%% �γ�Jacobi����
function FormJacobi %�γ��Ÿ��Ⱦ���,Hij����Pi��fj�󵼣�Nij����Pi��ej�󵼣�Jij����Qi��fj�󵼣�Lij����Qi��ej��
global Jacobi
global H N J L;
global BusNum;
global P Q;
global G B;
global Um Ua DeltaUa;
global PQnum PVnum;

mpc = case14;

%�����ʼ��
H = zeros(BusNum);
N = zeros(BusNum);
J = zeros(BusNum);
L = zeros(BusNum);
Jacobi = zeros(2*PQnum+PVnum);
inode = 2;
irow  = 1; % Jacobi�������
iclm  = 1; % Jacobi�������
%���ſ˱�ϵ��
for i=2:BusNum
    for j=2:BusNum
       if i==j
            H(i,i) = Q(i) + B(i,i) * Um(i) * Um(i);
            N(i,i) =-P(i) - G(i,i) * Um(i) * Um(i);
            J(i,i) =-P(i) + G(i,i) * Um(i) * Um(i);
            L(i,i) =-Q(i) - B(i,i) * Um(i) * Um(i); 
        else
            H(i,j) = -Um(i) * Um(j) * ( G(i,j) * sind(DeltaUa(i,j)) - B(i,j) * cosd(DeltaUa(i,j)) );
            N(i,j) = -Um(i) * Um(j) * ( G(i,j) * sind(DeltaUa(i,j)) + B(i,j) * cosd(DeltaUa(i,j)) );
            J(i,j) =  Um(i) * Um(j) * ( G(i,j) * cosd(DeltaUa(i,j)) + B(i,j) * sind(DeltaUa(i,j)) ); 
            L(i,j) = -Um(i) * Um(j) * ( G(i,j) * cosd(DeltaUa(i,j)) - B(i,j) * sind(DeltaUa(i,j)) );
        end    
    end 
end
%�γ��ſ˱Ⱦ���,�ⲿ��ע��Jacobi���󣬲ο���صļ�����ţ�����������㣬�����󵼵Ķ����Ϊ���ࡣPQ-PQ,PQ-PV,PV-PQ,PV-PV
row=1;
rowdelta=0;
for i=2:BusNum
    itype = mpc.bus(i,2);
    line=1;
    row = row + rowdelta;
    for j=2:BusNum
        jtype = mpc.bus(j,2); 
        if itype == 1
            rowdelta = 2;
            if jtype == 1   %PQ-PQ
                Jacobi(row,line)     = H(i,j);
                Jacobi(row,line+1)   = N(i,j);
                Jacobi(row+1,line)   = J(i,j);
                Jacobi(row+1,line+1) = L(i,j);
                line = line + 2;
            elseif jtype ==2 %PQ-PV
                Jacobi(row,line)     = H(i,j);
                Jacobi(row+1,line)   = J(i,j);
                line = line + 1;
            end
        elseif itype == 2
            rowdelta = 1;
            if jtype == 1 %PV-PQ
                Jacobi(row,line)     = H(i,j);
                Jacobi(row,line+1)   = N(i,j);
                line = line + 2;
            elseif jtype == 2 %PV-PV
                Jacobi(row,line)     = H(i,j);
                line = line + 1;
            end
        end
    end
end
end
%% ����������̵õ���ѹ����
function GetRevised
global Jacobi;
global deltaUm deltaUa;
global Um Ua;
global Delta; %�ڵ㹦�ʵĲ�ƽ����
global BusNum;
mpc = case14;
line=1;
deltaUm = zeros(BusNum,1);
deltaUa = zeros(BusNum,1);

deltaU = Jacobi \ Delta;
for i=2:BusNum
    itype = mpc.bus(i,2);
    if itype==1
        deltaUa(i) = deltaU(line);
        deltaUm(i) = deltaU(line+1) * Um(i);
        line = line + 2;
    elseif itype==2
        deltaUa(i) = deltaU(line);
        line = line + 1;
    end
end
end
%% �õ���ѹ��ֵ
function GetNewVoltageValue
global BusNum;
global deltaU deltaUm deltaUa;
global Um Ua;

for i = 2:BusNum
    Um(i) = Um(i) - deltaUm(i);
    Ua(i) = Ua(i) - deltaUa(i);
end

%���Ƕ�ֵ������0~360��
for i=1:BusNum
    if Ua(i)>=360 || Ua(i)<=-360
        A = Ua(i);
        B = fix(Ua(i));
        C = mod(B,360);
        Ua(i) = C + A - B;
    end
end
end
%% ����ƽ��ڵ㹦��
%% �ж��Ƿ񳬹�Լ������


%% ���

