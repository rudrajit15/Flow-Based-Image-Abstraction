function [image_op] = LineExtraction2(img,sigma_m,sigma_c)
%img = imread('4.2.04.tiff');
%img = rgb2gray(img);
sz = size(img);
sigma_s = 1.6*sigma_c;
%T = floor(10*sigma_c);
%S = floor(10*sigma_m);
S = 8;
T = 4;
max_ST = max(S,T);
win = 2*max_ST + 1;
winby2 = floor(win/2);
[t,image] = ETF2(img,win);
%image = imgaussfilt(image,0.5);
sz2 = size(image);
%[Tx, Ty, image] = myETF(img,win,0.2);
%sz_t = size(Tx);
%s_t = [sz_t(1) sz_t(2) 2];
%t = zeros(s_t);
%t(:,:,1) = Tx;
%t(:,:,2) = Ty;
%[image , t(:,:,1),t(:,:,2)] = myETF(img,11,1);
tau1 = 0.95;
tau2 = 0;
rho = 0.99;
g(:,:,1) = t(:,:,2);
g(:,:,2) = -t(:,:,1);

[gmag ,gdir] = imgradient(g(:,:,1),g(:,:,2));
[tmag ,tdir] = imgradient(t(:,:,1),t(:,:,2));
        
for run = 1:6
    H_g = zeros(size(image));
    for v= winby2+1: winby2+sz(2)
        for u = winby2+1:winby2+sz(1)
            angle = gdir(u,v);
            sum_t = 0;
                
            if ((angle >= -22.5) && (angle <= 22.5))
                for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u+q,v)*f_t;
                end
                
            elseif ((angle <=-157.5) || (angle >= 157.5)) 
                for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u-q,v)*f_t;
                end
                
            elseif ((angle >= 22.5) && (angle <= 67.5))
                for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u+q,v+q)*f_t;
                end
                
            elseif ((-157.5 <= angle) && (angle  <= -112.5))
                for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u-q,v-q)*f_t;
                end

            elseif ((angle >= 67.5) && (angle <= 112.5))
                 for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u,v+q)*f_t;
                 end
                 
            elseif ((angle <= -67.5) && (angle >= -112.5))
                 for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u,v-q)*f_t;
                end

            elseif ((112.5 <= angle) && (angle <= 157.5))
                 for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u-q,v+q)*f_t;
                 end
                 
            elseif ((angle >= -67.5) && (angle <= -22.5))
                for q = -T:T
                    f_t = normpdf(q,0,sigma_c) - rho*normpdf(q,0,sigma_s);
                    sum_t = sum_t + image(u+q,v-q)*f_t;
                 end
            end
            H_g(u,v) = sum_t;
        end
    end
    
    H_e = zeros(size(image));
    for v= 1:sz2(2)-2*max_ST-1
        for u = 1:sz2(1)-2*max_ST-1
            u2 = u;
            v2 = v;
            xc = u;
            yc = v;
            sum_s = 0;
            
            for f=1:2*S+1
                angle = tdir(u2,v2);
                
                if(f == S+1)
                    xc = u2;
                    yc = v2;
                end;
                
                g_s = normpdf(f-S-1,0,sigma_m);
                sum_s = sum_s + H_g(u2,v2)*g_s;
                
                if ((angle >= -22.5) && (angle <= 22.5)) 
                        %sum_s = sum_s + H_g(u2+1,v2)*g_s;
                        u2 = u2+1;
                        
                elseif ((angle <=-157.5) || (angle >= 157.5))
                        %sum_s = sum_s + H_g(u2-1,v2)*g_s;
                        u2 = u2-1;
                        
                elseif ((angle >= 22.5) && (angle <= 67.5))
                        %sum_s = sum_s + H_g(u2+1,v2+1)*g_s;
                        u2 = u2+1;
                        v2 = v2+1;

                elseif (-157.5 <= angle <= -112.5)
                        %sum_s = sum_s + H_g(u2-1,v2-1)*g_s;
                        u2 = u2-1;
                        v2 = v2-1;
                        
                elseif ((angle >= 67.5) && (angle <= 112.5))   
                        %sum_s = sum_s + H_g(u2,v2+1)*g_s;
                        v2 = v2+1;
                        
                elseif ((angle <= -67.5) && (angle >= -112.5))
                        %sum_s = sum_s + H_g(u2,v2-1)*g_s;
                        v2 = v2-1;
                        
                elseif (angle >= -67.5) && (angle <= -22.5)  
                        %sum_s = sum_s + H_g(u2+1,v2-1)*g_s;
                        u2 = u2+1;
                        v2 = v2-1;
                        
                elseif (112.5 <= angle && angle <= 157.5)
                        %sum_s = sum_s + H_g(u2-1,v2+1)*g_s;
                        u2 = u2-1;
                        v2 = v2+1;
                end;
            end
            
            H_e(xc,yc) = sum_s;
        end
    end
    
    linedimg = zeros(size(image));
    for v= winby2 +1: winby2 +sz(2)
        for u = winby2+1:winby2+sz(1)
            if((H_e(u,v) < 0 ) && (1+tanh(H_e(u,v)) < tau1))
            %if((H_e(u,v) < 0 ) || (tanh(H_e(u,v)) < tau2))
                linedimg(u,v) = 0;
            else
                linedimg(u,v) = 255;
            end
        end
    end
    
    image = (image.*linedimg)/255;
end
s = size(img);
image_op = zeros(s);
i = 1:s(1);
j = 1:s(2);
image_op(i,j) = linedimg(i+winby2,j+winby2);
%imshow(uint8(image_op));


