clear;
clc;

rng(144);

load gyyt_navi_data
load AKF_result.mat

% gyyt
gyyt_n = run_gyyt_filter([(kf_sol(:,1)+akf_sol(:,1))/2 kf_sol(:, 4)]');
gyyt_e = run_gyyt_filter([akf_sol(:,2) kf_sol(:, 5)]');

% sage-huga KF
[SH_N,~] = run_sage_huga_kf(kf_sol(:,[1,4])');
[SH_E,~] = run_sage_huga_kf(kf_sol(:,[2,5])');

% rtk
U0 = 1.0e+06 *[-2.179750330532283; 4.383964569544397; 4.074104603047358];
CLIGHT = 299792458.0; D2R = pi/180;
v0 = [0; 0; 0];

[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(U0,v0);
B0 = L_b; L0 = lambda_b;

T  = [ -sin(B0)*cos(L0) -sin(B0)*sin(L0) cos(B0)
    -sin(L0)          cos(L0)         0
    cos(B0)*cos(L0)  cos(B0)*sin(L0) sin(B0)];

panda_file = '20191219-rtk.dat';
p_llh = load_panda_llh(panda_file);
pd_local = zeros(size(p_llh,1),3);
for k = 1:size(p_llh,1)
    lat = p_llh(k,2)*pi/180;
    lon = p_llh(k,3)*pi/180;
    alt = p_llh(k,4);
    U_tmp = BLH2XYZ(lat, lon, alt);
    pd_local(k,:) = (T*(U_tmp-U0))';
end

figure;
hold on;
plot(pd_local(105:end, 1), pd_local(105:end,2), 'o-');
plot(lsq_sol(:,1), lsq_sol(:,2), '*-');
plot(kf_sol(:,1), kf_sol(:,2), 'd-');
plot(SH_N(1,:), SH_E(1,:), '.-');
plot(gyyt_n(1, :),  gyyt_e(1,:), '+-');
% xlabel('东向(m)');
% ylabel('北向(m)');
% legend('参考路径(RTK)','最小二乘', '卡尔曼滤波','Sage滤波','广义延拓滤波');
xlabel('East/m');
ylabel('North/m');
legend('Reference(RTK)','Lsq', 'Kalman Filter', 'Sage Filter', 'Generalized Extension Filter','Location', 'SouthEast');
print('-djpeg', '-r600', 'fig8.jpeg');
% 
xy0 = pd_local(105:105+40,:);
xy1=pd_local(105+41:105+149,:);
xy2=pd_local(105+150:105+310,:);
xy3 = pd_local(105+311:105+640,:);
xy4 = pd_local(105+641:end,:);

x0 = linspace(1,220, size(xy0,1));
fit_xy0(1,:) = spline(x0, xy0(:,1), 1:220);
fit_xy0(2,:) = spline(x0, xy0(:,2), 1:220);
x1 = linspace(1,110,size(xy1,1));
fit_xy1(1,:) = spline(x1,xy1(:,1),1:110);
fit_xy1(2,:) = spline(x1,xy1(:,2),1:110);
% figure; plot(x1, xy1(:,1), '.-', 1:110, fit_xy1(1,:), '*-'); hold on; plot(1:110, kf_sol(1:110,1),'d-');
% figure; plot(x1, xy1(:,2),'.-',  1:110, fit_xy1(2,:), '*-'); hold on; plot(1:110, kf_sol(1:110,2),'d-');
x2 = linspace(1,200,size(xy2,1));
fit_xy2(1,:) = spline(x2,xy2(:,1),1:200);
fit_xy2(2,:) = spline(x2,xy2(:,2),1:200);
x3 = linspace(1,320,size(xy3,1));
fit_xy3(1,:) = spline(x3,xy3(:,1),1:320);
fit_xy3(2,:) = spline(x3,xy3(:,2),1:320);
x4 = linspace(1,1118-850,size(xy4,1));
fit_xy4(1,:) = spline(x4,xy4(:,1),1:(1118-850));
fit_xy4(2,:) = spline(x4,xy4(:,2),1:(1118-850));

pd_fit =  [fit_xy0 fit_xy1 fit_xy2 fit_xy3 fit_xy4]';

% figure; plot(pd_fit(:,1), '.-'); hold on; plot(kf_sol(:,1),'o-'); legend('fit', 'KF x');
% figure; plot(pd_fit(:,2), '.-'); hold on; plot(kf_sol(:,2),'o-'); legend('fit', 'KF y');

Nx = 1100;
rms_lsq = norm(sqrt((pd_fit(1:Nx,:) - lsq_sol(1:Nx,1:2)).^2/Nx))
rms_kf = norm(sqrt((pd_fit(1:Nx,:) - kf_sol(1:Nx,1:2)).^2/Nx))
rms_gyyt = norm(sqrt((pd_fit(1:Nx,:) - [gyyt_n(1,1:Nx); gyyt_e(1,1:Nx)]').^2/Nx))
rms_sage = norm(sqrt((pd_fit(1:Nx,:) - [SH_N(1,1:Nx); SH_E(1,1:Nx)]').^2/Nx))
rms_lsq2 = sqrt(sum((pd_fit(1:Nx,1) - lsq_sol(1:Nx,1)).^2+(pd_fit(1:Nx,2) - lsq_sol(1:Nx,2)).^2))/sqrt(Nx)

% 生成

function gx = run_gyyt_filter(x)

  dt = 1.0;
  
  P = diag([2 0.2])*5;
  R = diag([1.^2 0.1^2]);
  H = [1 0; 0 1];
  q = 0.8;
  F = [0 1;
       0 0];
  [A,Q] = lti_disc(F,[],diag([0 q]),dt);
  N = 15;
  
  gx = x(:,1:N);
  for k=N+1:size(x,2)
      
      if k == N+1
          % Init the gyyt filter
          
          prev_x = x(:,k-N:k-1);
          prev_P = P;
          gyyt_filter = GyytFilter3(prev_x, prev_P, A, H, Q, R, N,2,0);
      end
      
      if k >= N+1
          % perform the filter
          gyyt_filter.update(x(:,k));
          
          gyyt_filter.P0 = gyyt_filter.est_P;
          
          gyyt_filter.prev_x = [gyyt_filter.prev_x(:, 2:end) gyyt_filter.est_x];
          
          gx(:,k) = gyyt_filter.est_x;
          
      end
  end
end

function [X, P] = run_sage_huga_kf(x)
  dt = 1.0;
  P = diag([2 0.2]);
  R = diag([1.^2 0.1^2])*0.1;
  H = [1 0; 0 1];
  q = 0.5;
  F = [0 1;
       0 0];
  [A,Q] = lti_disc(F,[],diag([0 q]),dt);
  G = eye(size(Q));
  
  X0 = x(:,1);
  Z = x;
  
  [X,P]=Sage_HusaKF(A,G,H,Q,R,X0,Z,P)
   
end