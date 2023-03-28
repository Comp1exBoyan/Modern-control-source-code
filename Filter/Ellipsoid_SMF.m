% The application of Set Membership Filter in localizattion 
% Author: Yifeng Tang ^ University of Chinese Academy of Sciences
% Introduction: a kind of SMF based on ellipsoid
function Ellipsoid_SMF
clc;clear;
clc;clear;
T=1;            % sample peroid
N=30/T;         % sample times 
X=zeros(4,N);   % real position and velocity 
X(:,1)=[-100,2,200,20];%initial pos and v 
Z=zeros(2,N);   % measurement value
Z(:,1)=[X(1,1),X(3,1)];
delta_w=1e-2;
Q=delta_w*diag([0.5,1,0.5,1]);
R=100*eye(2);% mean value of system noise
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];% transformation matrix
H=[1,0,0,0;0,0,1,0];% measurement matrix
G =eye(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:N
    X(:,t)=F*X(:,t-1)+sqrtm(Q)*randn(4,1);% the real trajectory of the obj 
    Z(:,t)=H*X(:,t)+sqrtm(R)*randn(2,1); %  measurement 
end
% Set membership 
Pk = eye(4);
Xkf=zeros(4,N);
Xkf(:,1)=X(:,1);% initialization
for i=2:N
        syms pk;
        syms qk ;
    % time updating
    Xn = F*Xkf(:,i-1);
    eig_mat = G*Q*G'*((F*Pk*F')^(-1));
    eig_p = eig(eig_mat);
    eq_pk = 1/(eig_p(1)+pk) + 1/(eig_p(2)+pk) + 1/(eig_p(3)+pk) + 1/(eig_p(4)+pk)  == 4/(pk^2+pk);
    pk = solve(eq_pk);
    pk = double(pk(1));
    P_trans = subs((pk+1)*F*Pk*F' + (1/pk + 1)*G*Q*G');
    ek = Z(:,i) - H* Xn;
    
    % measurament updating
    Sk = qk*H*P_trans*H' + R;
    beta = 1 + qk-qk*ek'*Sk^(-1)*ek;
    eig_matt = P_trans *H'*R^(-1)*H;
    eig_q = eig(eig_matt);
    eq_qk = eig_q(1)/(1+eig_q(1)*qk) +  eig_q(2)/(1+eig_q(2)*qk) +eig_q(3)/(1+eig_q(3)*qk) + eig_q(4)/(1+eig_q(4)*qk) == 4 * 1/beta *(1 - ek'*Sk^(-1)*R*Sk^(-1)*ek);
    qk = vpa(solve(eq_qk));
    qk = max(qk);
    if 4*(1-ek'*R^(-1)*ek) - trace(P_trans*H'*R^(-1)*H) >=0
        qk = 0;
    end
    
    % filtering updating
    Xkf(:,i) = Xn+ subs(qk*P_trans*H'*Sk^(-1)*ek);
    Pk = subs(beta)*(P_trans - subs(qk*P_trans * H' *Sk^(-1)*H*P_trans));
end

for i=1:N
    Err_Observation(i)=RMS(X(:,i),Z(:,i));% errors before filtering 
    Err_KalmanFilter(i)=RMS(X(:,i),Xkf(:,i));% after filtering 
end
%画图
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k');%real trajectory
plot(Z(1,:),Z(2,:),'-b.');%measurement trajectory
plot(Xkf(1,:),Xkf(3,:),'-r+');%filtered trajectory
legend('real trajectory','measurement trajectory','filtered trajectory')
figure
hold on; box on;
plot(Err_Observation,'-ko','MarkerFace','g')
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
legend('before','after')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%norm 2 -- Euclidean distance
function dist=RMS(X1,X2)
if length(X2)<=2
    dist=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );
else
    dist=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%
