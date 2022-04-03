clc
clear
fluid1='oxygen';
refpropm('T','C',0,'',0,fluid1);
refpropm('P','C',0,'',0,fluid1);
refpropm('Z','C',0,'',0,'oxygen')
syms x z
solve(x-5==0,x)
y=0.00001;
for i=1:10000
    e=1/(y-0.07780)-0.45724/y/(1.07780*y+0.0717)-1;
    y=y+0.0001;
end
solve(1/(x-0.07780)-0.45724/x/(1.07780*x+0.0717)-1==0)
