classdef GyytFilter3 < handle
   
    properties
       prev_x;
       P0;
       A;
       H;
       
       N;
       
       n;
       
       pred_x;
       est_x;
       
       pred_P;
       est_P;
       
       Q;
       R;
       
       debug;
       kf_pred;
    end
    
    methods
        function obj = GyytFilter3(prev_x_, prev_P_, A_, H_, Q_, R_, N_, n_, debug_)
            
            if nargin < 8
                n_ = 2;
            end
            
            obj.prev_x = prev_x_;
            obj.P0 = prev_P_;
            obj.H = H_;
            obj.Q = Q_;
            obj.R = R_;            
            obj.A = A_;
            obj.N = N_;
            
            obj.n = n_;
            
            obj.debug = debug_;
            obj.kf_pred = 0;
            
            obj.pred_x = [];
            obj.pred_P = [];
            obj.est_x = [];
            obj.est_P = [];
            
        end
       
        function gyyt_predict(obj)
            % get the interpolation coef
            coef = obj.gyyt_interp()';
            
            t = obj.N+1;
            
            % calculate the interpolation value
            if obj.n == 2
                obj.pred_x(1,:) = coef(1) + coef(2)*t + coef(3)*t.*t;
                obj.pred_x(2,:) = coef(2)+coef(3)*2*t;
            else
                obj.pred_x(1,:) = coef(1) + coef(2)*t + coef(3)*t.*t + coef(4)*t.^3;
                obj.pred_x(2,:) = coef(2)*t + coef(3)*2*t + coef(4)*3*t.^2;
            end
            
%             x0 = obj.prev_x(:,end-1);
%             x1 = obj.pred_x;
%             obj.A = x1 * x0';

            Pt = obj.P0;
            P = Pt(:,:,end); 
            
            P = 1.05^2*obj.A * P * obj.A' + obj.Q;
            
            obj.pred_P = P;

        end
        
        function coef = gyyt_interp(obj)
            
            coef = zeros(obj.n+1,1);
            
            x = obj.prev_x;
            
%             Px = obj.P0(1,1,:); Px(1:end-5)=10^8;
%             
%             [~, n_ind] = min(Px);
            
                        
            n_ind = obj.N;
            
            ti = 1:obj.N; xi = x(1, 1:obj.N); xdi = x(2, 1:obj.N);
            tn = n_ind;
            xn = x(1, n_ind);
            xdn = x(2, n_ind);
            ti(n_ind) = []; xi(n_ind) = []; xdi(n_ind) = [];
         
%             c11 = sum(ti.^2-2*ti*tn+tn.^2);
%             c12 = sum(ti.^3-ti.^2*tn-ti.*tn.^2+tn.^3);
%             c21 = c12;
%             c22 = sum(ti.^4-2*ti.^2*tn.^2+tn.^3);
            if obj.n == 2            
%                 c11 = sum(ti - tn);
%                 c12 = sum(ti.^2 - tn.^2);
%                 c21 = c12;
%                 c22 = sum((ti-tn).*((ti+tn).^2+4));

                c11 = sum(ti - tn);
                c12 = sum(ti.^2 - tn.^2);
                c21 = c12;
                c22 = sum((ti+tn).^2.*(ti-tn));

                vflag = 0;

                b1 = sum(xi-xn);
                b2 = sum((ti+tn).*(xi-xn)+2*(xdi-xdn)*vflag);

                C = [c11 c12; c21 c22];
                B = [b1; b2];
                cc = C \ B;

%                 alpha = 1;
%                 
%                 Ac = sum(2*alpha.*(xn-xi+(ti-tn).*xdn).*(ti-tn)+4.*(xdn-xdi));
%                 Bc = sum(2*alpha*(ti-tn).^3+8*(ti-tn));
%                 
%                 a2 = -Ac/Bc;
%                 
%                 a1 = xdn - 2*a2*tn;

                a1 = cc(1); a2 = cc(2);
                
                a0 = xn - a1*tn - a2*tn^2;

                coef = [a0; a1; a2];

            else
                
                c11 = sum(ti.^2-2.*ti.*tn+tn.^2);
                c12 = sum(ti.^3-ti.^2.*tn-tn.^2.*ti+tn.^3);
                c13 = sum(ti.^4-ti.^3.*tn-tn.^3.*ti+tn.^4);
                c21 = c12;
                c22 = sum(ti.^4-2*ti.^2*tn.^2+tn.^4);
                c23 = sum(ti.^5-ti.^2*tn.^3-tn.^2*ti.^3+tn.^5);
                c31 = c13;
                c32 = c23;
                c33 = sum(ti.^6-2*ti.^3*tn.^3+tn.^6);
                
                b1 = sum(ti.*xi-ti.*xn-tn.*xi+tn.*xn);
                b2 = sum(ti.^2.*xi-ti.^2.*xn-tn.^2.*xi+tn.^2.*xn);
                b3 = sum(ti.^3.*xi-ti.^3*xn-tn.^3.*xi+tn.^3.*xn);
                
                C = [c11 c12 c13; c21 c22 c23; c31 c32 c33];
                B = [b1; b2; b3];
                
                cc = C\B;
                
                a = xn - cc(1)*tn - cc(2)*tn^2 - cc(3)*tn^3;
                
                coef = [a; cc(1); cc(2); cc(3)];
            end
            
%             coef1 = polyfit(1:obj.N, x(1,:), obj.n);
%             
%             coef = coef1(end:-1:1);
%             
            if obj.debug
               figure(100); hold on;
               plot(x(1,:),'o-'); hold on;
               tt = 1:0.1:obj.N+1;
               if obj.n == 2
                    yy = coef(1) + coef(2)*tt + coef(3)*tt.*tt;
               else
                   yy = coef(1) + coef(2)*tt + coef(3)*tt.^2 + coef(4)*tt.^3;
               end
               plot(tt, yy, '*-');
               plot(tt(end), yy(end), 'bd');
               plot(tt(end), obj.kf_pred, 'kd');
               legend('用于外推的值','最小二乘拟合','最小二乘外推点','KF 预测值');
               obj.debug = obj.debug + 1;
               
               if obj.debug > 15
                   obj.debug = 0;
               end
               
               figure(1);
            end
        end
              
        
        function [X,P,K,IM,IS,LH] = gyyt_update(obj, X,P,y,H,R)
            
            %
            % Check which arguments are there
            %
            if nargin < 5
                error('Too few arguments');
            end
            
            %
            % update step
            %
            IM = H*X;
            IS = (R + H*P*H');
            K = P*H'/IS;
            X = X + K * (y-IM);
            P = P - K*IS*K';
            if nargout > 5
                LH = gauss_pdf(y,IM,IS);
            end
            
        end
         
        
        function update(obj, y)
            obj.gyyt_predict();
            [X,P] = obj.gyyt_update(obj.pred_x, obj.pred_P,y,obj.H,obj.R);
            obj.est_x = X;
            obj.est_P = P;
        end
    end
end

