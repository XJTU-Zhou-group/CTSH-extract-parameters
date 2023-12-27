clc
clear
p = gcp('nocreate');
if isempty(p)
    disp('启动并行运算，核心数：8');
    parpool('local', 8);
else
    disp(['并行运算已启动，核心数：' num2str(p.NumWorkers)]);
end
n=5;
low=[80080 6670 0.228 2930 2500];
high=[181000 10300 0.34 7170 3708];
% 设置非默认求解器选项
options = optimoptions("ga",'UseParallel', true,"Display","iter","PopulationSize",2000,"MaxGenerations",800,'PlotFcn', @gaplotbestf);
% 求解
[solution,objectiveValue] = ga(@errofdata,n,[],[],[],[],low,high,[],[],options);
sol(1,:)=solution;
writematrix(sol,"sol.csv")
writematrix(errofsol,"errofdata.csv");
str=num2str(1);
saveas(gcf,str,'fig'); %保存当前窗口的图像



