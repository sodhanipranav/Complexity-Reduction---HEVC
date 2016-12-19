clc
% ----- Run Once -----
 
clear store
%figure, imshow(y1{4}/255)
for k=1:2:n
    clear val
    clear claim
    clear nullset
    clear bmat
    clear big
    clear b_avg
    clear th
    y_best(:,:) = y1{k};
    y_rough(:,:) = y1{k+1};
    f_energy = 0; % frame energy
    b_energy = 0; % block energy
    count  = 0;
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    count5 = 0;
    count9 = 0;
    
    % calculate energy in the residual frame
    for i=1:64*16
        for j=1:64*30
            f_energy = f_energy + y_rough(i,j)^2;
        end
    end
    
    % using y_best, find out CTUs which have 64x64 size as the best mode
    % val contains CTU number of all blocks chosen
    for numCTU=1:480
        row = floor((numCTU-1)/30);
        col = rem((numCTU-1),30);
        
        if y_best((row+1)*64-32, (col+1)*64-32)~=255
            count=count+1;
            val(count)=numCTU;
        end
    end
    if count>0
    [temp, a]=size(val);
    best(k)=a;
    % calculate block energy of all 32 x 32 blocks using the residual obtained
    % as a result of rough motion estimation (use y_rough: devoid of boundaries)
    
    for numCTU=1:480*4
        b_energy=0;
        row = floor((numCTU-1)/60);
        col = rem((numCTU-1),60);
        
        for i=1+row*32:32+row*32
            for j=1+col*32:32+col*32
                b_energy = b_energy + y_rough(i,j)^2;
            end
        end
        b_avg(numCTU) = b_energy/(32*32);
    end
    
    x=1;
    y=60;
    
    while(y<480*4)
        
        for numCTU=x:2:y
            count2 = count2 + 1;
            big(count2) = b_avg(numCTU)+b_avg(numCTU+1) + b_avg(numCTU+60)+ b_avg(numCTU+61);
            
        end
        x=x+120;
        y=y+120;
    end
        
    % creating bmat breaking 64x64 into 32x32 blocks
    % contains residual energy info of each 32x32 block
    m = 60;
    temp=1;
    for numCTU=1:480*4
        if numCTU==m
            bmat(temp, rem(numCTU-(temp-1)*60, 61))=b_avg(numCTU);
            temp=temp+1;
            m=m+60;
        else
            bmat(temp, rem(numCTU-(temp-1)*60, 61))=b_avg(numCTU);
        end
    end
    count9=0;
    for i=1:480
        if big(i)<=median(big(val))
            count9 = count9 + 1;
            th(count9)=i;
        end
    end
    [temporary, len] = size(b_avg);
    %     uu=0;
    %     clear cost
    %     for thresh=3:0.1:13
    %         uu=uu+1;
    %         clear claim
    %         clear nullset
    %         count3 = 0;
    %         count5=0;
    
    for i=1:32 % 14
        for j=1:60 %26
            numCTU = (floor((i-1)/2))*30 + floor((j+1)/2);
            if big(numCTU)>mean(big) || bmat(i,j)>(median(b_avg))^1.2
               count3 = count3+1;
               nullset(count3) = numCTU;
            end
        end
    end
       
    claim =1:480;
    [temp, c] = size(nullset);
    for i=1:480
        for j=1:c
            if claim(i)==nullset(j)
                claim(i)=0;
            end
        end
    end
    
    % claim ready for comparison with val
    match=0;
    for i=1:a
        for j=1:480
            if val(i)==claim(j)
                match = match+1;
            end
        end
    end
    
    for i=1:480
        if claim(i)==0
            count5=count5+1;
        end
    end
    
    fa(k)=480-count5-match;
    fr(k)=a-match;
    
    my(k)=match;
    fenergy(k)=f_energy/100000;
    medbavg(k)=median(b_avg);
    less(k)=count9;
    end
end

clear store
store(1,:)=fenergy(:);
store(2,:)=best(:);
store(3,:)=my(:);
store(4,:)=fa(:);
store(5,:)=fr(:);
store(6,:)=medbavg(:);
store(7,:)=less(:);
sum(best)
sum(my)
sum(fa)
sum(fr)

%         cost(uu)=0.55*fa(k) + 0.45*fr(k);
% end
%     vec=3:0.1:13;
%     figure, plot(vec, cost);

% ----------------------

% for numCTU=1:91*4*4
%     b1_energy=0;
%     row = floor((numCTU-1)/52);
%     col = rem((numCTU-1),52);
%     
%     for i=1+row*16:16+row*16
%         for j=1+col*16:16+col*16
%             b1_energy = b1_energy + y_rough(i,j)^2;
%         end
%     end
%     b1_avg(numCTU) = b1_energy/(32*32);
%     
% end
% e = 52;
% temp1=1;
% for numCTU=1:91*4*4
%     if numCTU==e
%         b1mat(temp1, rem(numCTU-(temp1-1)*52, 53))=b1_avg(numCTU);
%         temp1=temp1+1;
%         e=e+52;
%     else
%         b1mat(temp1, rem(numCTU-(temp1-1)*52, 53))=b1_avg(numCTU);
%     end
% end
% 
% for i=1:28
%     for j=1:52
%         numCTU = (floor((i-1)/4))*13 + floor((j+3)/4);
%         if big(numCTU)>mean(big) || b1mat(i,j)>(median(b1_avg(th))) || bmat(floor((i+1)/2),floor((j+1)/2))>(mean(b_avg)) 
%             count3 = count3+1;
%             nullset(count3) = numCTU;
%         end
%     end
% end

% % ----------------------
%
% load rd.dat
% total=0;
%  for i=1:91
%      total = total + rd(i);
%  end
%  for numCTU=1:91
%      rd(numCTU) = (rd(numCTU)/total)*100;
%  end
%


