qpskConstellation = [-1+1i 1+1i; -1-1i 1-1i]/sqrt(2);
qpsk = reshape(qpskConstellation,1,[]); 
Num  = 40;
outter = 60;
 for nn = 1:outter
  qpsk = qpsk * (outter-1)/outter;
  c = rand(Num,3);       %随机生成了12种颜色。RGB随机。
     for idx = 1:Num
         theta = pi/N/Num*idx;
         rou = [cos(theta) sin(theta);sin(theta) -cos(theta)];
     realPart = real(qpsk);
     imagPart = imag(qpsk);
     reim = rou * [realPart;imagPart];
     realPart2 = real(qpsk*0.3);
     imagPart2 = imag(qpsk*0.3);
     reim2 = rou * [realPart2;imagPart2]; 
     plot(reim(1,:),reim(2,:),'o','color',c(idx,:))
     hold on;
     plot(reim2(1,:),reim2(2,:),'.','color',c(idx,:))
     hold on;
     pause(0.005);
     end
 end