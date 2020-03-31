% GuangYiYanTuo-based Kalman Filter demonstration with sine signal.
%
% History:
%    3.12.2002 SS  The first implementation
%
% Copyright (C) 2002 Simo S鋜kk?
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
clear;
clc;

set_color_noise = 1;

rng(144);

MT = 1000;
for mt = 1:MT
  %
  % Create sine function
  %
  S1 = [0.2;1.0];
  S2 = [1.0;-0.2];
  sd = 0.1;
  dt = 0.1;
  w = 1;
  T = (0:dt:30);
  X = sin(w*T);
  
  L = size(X, 2);
  c = [1 -0.9];
  nc = length(c)-1;
  xik = zeros(nc, 1);
  xi = randn(L, 1);
  
  e = zeros(L, 1);
  for k = 1:L
    e(k) = c*[xi(k); xik];
    for i=nc:-1:2
       xik(i) = xik(i-1); 
    end
    xik(1) = xi(k);
  end
  
%   X = (w*T).^2 + 2*(w*T) + 1;

  % 添加有色噪声
  if set_color_noise == 1
        Y = X + sd*randn(size(X)) + 0.1*e';
  else
  % 无有色噪声
        Y = X + sd*randn(size(X));
  end
  

  %
  % Initialize KF to values
  %
  %   x = 0
  %   dx/dt = 0
  %
  % with great uncertainty in derivative
  %
  M = [0;0];
  P = diag([0.1 2]);
  R = sd^2;
  H = [1 0];
  q = 0.1;
  F = [0 1;
       0 0];
  [A,Q] = lti_disc(F,[],diag([0 q]),dt);

    %
  % Track and animate
  %
  MM = zeros(size(M,1),size(Y,2));               GMM=MM; 
  GMP = MM; KMP = MM;
  PP = zeros(size(M,1),size(M,1),size(Y,2));     GPP=PP;
  clf;
  clc;
  
  N = 20;
  
  for k=1:size(Y,2)
      %
      % Track with KF
      %
            
      [M,P] = kf_predict(M,P,A,Q);
      x_pred = M(1,1);
      KMP(:,k) = M;
      [M,P] = kf_update(M,P,Y(k),H,R);
      
      MM(:,k) = M;
      PP(:,:,k) = P;
      
      if k == N+1
          % Init the gyyt filter
         GMM = MM;
         GMP(:,1:k) = KMP(:,1:k);
         prev_x = MM(:,k-N:k-1);
         prev_P = PP(:,:,k-N:k-1); GPP(:,:,1:N) = PP(:,:,k-N:k-1);
         gyyt_filter = GyytFilter3(prev_x, prev_P, A, H, Q, R, N,2,0); 
         
      end
      
      gyyt_filter.kf_pred = x_pred;
      
      if k >= N+1
          % perform the filter
          gyyt_filter.update(Y(k));
          GMM(:,k) = gyyt_filter.est_x;
          GPP(:,:,k) = gyyt_filter.est_P;
          
          gyyt_filter.P0 = GPP(:,:,k-N+1:k);
          
          gyyt_filter.prev_x = GMM(:,k-N+1:k);
          
          GMP(:,k) = gyyt_filter.pred_x;
          
          if k == N+1 && gyyt_filter.debug
             figure(100);
             hold on;
             plot(N+1, M(1,1), 'ro');
             plot(N+1, GMM(1,k), 'bo-');
             legend('用于外推的值','最小二乘拟合','最小二乘外推点','KF 预测值','KF更新值','最小二乘更新值');
             figure(1);
          end
      end
      
  end
  plot(T, X, 'o-'); hold on;
  plot(T, Y, 'o'); hold on;
  plot(T, MM(1,:), '*-');
  plot(T, GMM(1,:), 'd-');
%   plot(T, KMP(1,:), 'd-');
%   plot(T, GMP(1,:), 'd-');
  % legend('Real data','Sim data', 'KF', '广义延拓');
  legend('Real data', 'Simulation data', 'Kalman Filter', 'Generalized Extension Filter');
  xlabel('Time/s');
  ylabel('Y');
  
  if mt == 1
     if set_color_noise == 0
         print('-djpeg', '-r600', 'fig4.jpeg');
     else
         print('-djpeg', '-r600', 'fig6.jpeg');
     end
  end
  
 data = [T' X' Y' GMM(1,:)' MM(1,:)'];
 save 'data.txt' data '-ascii'
 
  KRMS = sqrt(mean((X-MM(1,:)).^2));
  GRMS = sqrt(mean((X-GMM(1,:)).^2));
  
  fprintf('RMS error:\n');
  fprintf('KF = %f\n', KRMS);
  fprintf('GMM = %f\n', GRMS);
  
  MRMS(mt, 1) = KRMS;
  MRMS(mt, 2) = GRMS;
  
end

figure; histogram(MRMS(:,1),100);
hold on;
% figure; 
histogram(MRMS(:,2),100);
% legend('KF','广义延拓');
legend('Kalman Filer', 'Generalized Extension Filter');
title(sprintf('L=%d',N));
xlabel('Error');
ylabel('Counts per bin')
if set_color_noise == 0
    print('-djpeg', '-r600', 'fig5.jpeg');
else
    print('-djpeg', '-r600', 'fig7.jpeg');
end