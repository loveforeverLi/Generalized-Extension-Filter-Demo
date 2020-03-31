function [X,P]=Sage_HusaKF(F,G,H,Q,R,X0,Z,P)
% Sage-Husa adeptive KF 
Y = []; 
N=length(Z); 
M=length(X0); 
X=zeros(M,N); 
X(:,1)=X0; 
s=1*eye(2); 
q= zeros(M,1); 
r = 0;  
b = 0.970;
for k=2:N 
    X_est=F*X(:,k-1)+q;                  %����һ��Ԥ����ƣ�X(k/k-1) 
    P_pre=F*P*F'+G*Q*G';               %һ��Ԥ����Ƶľ������P(k/k-1) 
    e(:,k)=Z(:,k)-H*X_est-r;           %����в�epsilon(k) 
    K=P_pre*H'*inv((H*P_pre*H')+R);    %kʱ�̵������� 
    X(:,k)=X_est+K*e(:,k);           %kʱ�̵�״̬����X(k) 
    P = (eye(M)-K*H)*P_pre*(eye(M)-K*H)'+K*R*K';  %����������P(k) 
 
%       r = 1/k*((k-1)*r +Z(:,k)-H*X_est);     
%       q = 1/k*((k-1)*q+X(:,k)-F*X(:,k-1));      
%       R = 1/k*((k-1)*R+Z(:,k)*Z(:,k)'-H*P*H');      
%       Q = 1/k*((k-1)*Q+K*e(:,k)*e(:,k)'*K'+P-F*P*F'); 
   d = (1-b)/(1-b^(k)); 
   % v = Z(:, k)-H*X_est;
    r = (1-d)*r +d*(Z(:,k)-H*X_est); 
    q = (1-d)*q +d*(X(:,k)-F*X(:,k-1)); 
    R = (1-d)*R +d*(e*e'-H*P*H'); 
    % Q = (1-d)*Q +d*(K*e(:,k)*e(:,k)'*K'+P-F*P*F');
     Y=[Y,Q];
end