%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Incidence
%data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',11,'C3:AG57');%55 cities within 200km to the center cities
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');%200 多 cities
%w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',12);%55 cities within 200km to the center cities
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);%200 多 cities 
W = sparse(w);
nobs = length(data);
%%% for Incidence
y = data(:,1);%incidence(例/万人)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%sum + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
x = [xx x28 x30]; %sum_Ratio09 (normalized)
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)]; %sum_Ratio09
%x = [ones(nobs,1) data(:,6:18) data(:,29:30)]; %sum_Rank09
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Rank09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明可以拒绝原假设，则存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明不可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明不可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明不可以拒绝原假设
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实不显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。表明可用sac模型
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('incidence','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho变得不显著了
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho依然不显著

%综合以上情况，在单独的SEM中lambda不显著，但在加入spatial lag term of dependent
%variable（实质是变为类似于SAC模型）之后,lambda显得显著了。在单独的FAR中rho显著(rho=0.107994, p-value=0.000449)，但在加入explanatory
%variable之后，rho变得不显著了。因此在SAC中，lambda应是显著的，而rho应是不显著的。res2 = sac(y,x,W,W2);可以满足这一要求，虽然它的R方不是最大的。
%%SAC
W2 = slag(W,2); % standardized W2 result from slag(The function slag can generate a second order spatial weight matrix. A function slag can be used to produce higher order spatial lags.)
res1 = sac(y,x,W2,W);% general spatial model W2,W
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2 (selected)
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2
prt(res4,vnames); % print the output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Incidence
%The second round of variable selection

%data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',11,'C3:AG57');%55 cities within 200km to the center cities
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');%200 多 cities
%w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',12);%55 cities within 200km to the center cities
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);%200 多 cities 
W = sparse(w);
nobs = length(data);
%%% for Incidence
y = data(:,1);%incidence(例/万人)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%sum + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));%sum_Ratio09 (normalized)
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
%%x = [xx x28 x30]; 
x = [xx(:,1) xx(:,3:6) xx(:,12) x28 x30];%第一轮中回归系数显著的变量有Density,PGDP,Income,Hospital,CollegeStu,sum_Ratio,mc_distance
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)]; %sum_Ratio09
%x = [ones(nobs,1) data(:,6:18) data(:,29:30)]; %sum_Rank09
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
%%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
vnames = strvcat('Incidence','const','Density','PGDP','Income','Hospital','CollegeStu','sum_Ratio09','mc-distance');
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Rank09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明可以拒绝原假设，则存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明不可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明不可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明不可以拒绝原假设
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实不显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。表明可用sac模型
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('incidence','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
%%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
vnames = strvcat('Incidence','const','Density','PGDP','Income','Hospital','CollegeStu','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho变得不显著了
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho依然不显著

%综合以上情况，在单独的SEM中lambda不显著，但在加入spatial lag term of dependent
%variable（实质是变为类似于SAC模型）之后,lambda显得显著了。在单独的FAR中rho显著(rho=0.107994, p-value=0.000449)，但在加入explanatory
%variable之后，rho变得不显著了。因此在SAC中，lambda应是显著的，而rho应是不显著的。res2 = sac(y,x,W,W2);可以满足这一要求，虽然它的R方不是最大的。
%%SAC
W2 = slag(W,2); % standardized W2 result from slag(The function slag can generate a second order spatial weight matrix. A function slag can be used to produce higher order spatial lags.)
res1 = sac(y,x,W2,W);% general spatial model W2,W
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2 (selected)
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2
prt(res4,vnames); % print the output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Startweek (without closeness)
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);
W = sparse(w);
nobs = length(data);
%%% for Startweek
y = data(:,3);%startweek(09年的第几周)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%% closeness + sumR + mcd
%x = [ones(nobs,1) data(:,6:18) data(:,26) data(:,28) data(:,30)];
%%% sumR + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
x = [xx x28 x30]; %sum_Ratio09 (normalized)
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)];
%vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
%vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明可以拒绝原假设，则存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明可以拒绝原假设，lambda是显著的
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('Startweek','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho显著
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho不显著
%%SAC
W2 = slag(W,2); % standardized W2 result from slag
res1 = sac(y,x,W2,W);% general spatial model W2,W (selected)
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2
prt(res4,vnames); % print the output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Duration (without closeness)
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);
W = sparse(w);
nobs = length(data);
%%% for Duration
y = data(:,4);%duration(周数)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%% closeness + sumR + mcd
%x = [ones(nobs,1) data(:,6:18) data(:,26) data(:,28:29)];
%%% sumR + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
x = [xx x28 x30]; %sum_Ratio09 (normalized)
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)];
%vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
%vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明可以拒绝原假设，则存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明可以拒绝原假设
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('Duration','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho变得不显著了
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho依然不显著
%%SAC
W2 = slag(W,2); % standardized W2 result from slag
res1 = sac(y,x,W2,W);% general spatial model W2,W
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2 (selected)
prt(res4,vnames); % print the output 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Startweek (with closeness)
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);
W = sparse(w);
nobs = length(data);
%%% for Startweek
y = data(:,3);%startweek(09年的第几周)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%% closeness + sumR + mcd
%x = [ones(nobs,1) data(:,6:18) data(:,26) data(:,28) data(:,30)];
%%% sumR + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x26=data(:,26);
x26=(x26-min(x26))./(max(x26)-min(x26));
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
x = [xx x26 x28 x30]; %sum_Ratio09 (normalized)
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)];
%vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
%vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','Closeness','sum_Ratio09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明不可以拒绝原假设，则不存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明不可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明不可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明不可以拒绝原假设，lambda是显著的
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('Startweek','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Startweek','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','Closeness','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho显著
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho不显著
%%SAC
W2 = slag(W,2); % standardized W2 result from slag
res1 = sac(y,x,W2,W);% general spatial model W2,W (selected)
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2
prt(res4,vnames); % print the output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Duration (with closeness)
data = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',10,'C3:AH275');
w = xlsread('E:\科研\Spatial statistics\input-output flows\econometrics\spatial econometrics manual\H1N1',9);
W = sparse(w);
nobs = length(data);
%%% for Duration
y = data(:,4);%duration(周数)
y = (y-min(y))./(max(y)-min(y)); %normalized
%%% closeness + sumR + mcd
%x = [ones(nobs,1) data(:,6:18) data(:,26) data(:,28:29)];
%%% sumR + mcd
xx=ones(nobs,1);
for i=6:18
    x=data(:,i);
    x=(x-min(x))./(max(x)-min(x));
    xx=[xx x];
end
x26=data(:,26);
x26=(x26-min(x26))./(max(x26)-min(x26));
x28=data(:,28);
x28=(x28-min(x28))./(max(x28)-min(x28));
x30=data(:,30);
x30=(x30-min(x30))./(max(x30)-min(x30));
x = [xx x26 x28 x30]; %sum_Ratio09 (normalized)
%x = [ones(nobs,1) data(:,6:18) data(:,28) data(:,30)];
%vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','road cap','railway cap','flight cap','degree_1','degree_2','degree_3','betweenness','closeness','EigenCen','sum_Ratio09','mc-distance');
%vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','Closeness','sum_Ratio09','mc-distance');
spy(W);
%检验spatial correlation in the least-squares residuals of the regression model是否存在
%如果检验结果表明spatial correlation in the least-squares residuals存在，那么the SEM model would be an appropriate way to proceed
%零假设是no spatial correlation exists in the least-squares residuals of the regression model
res1 = moran(y,x,W);
prt(res1);%结果表明可以拒绝原假设，则存在spatial correlation in the least-squares residuals of the regression model，适用SEM model
res2 = lmerror(y,x,W);
prt(res2);%结果表明不可以拒绝原假设
res3 = lratios(y,x,W);
prt(res3);%结果表明不可以拒绝原假设
res4 = walds(y,x,W);
prt(res4);%结果表明可以拒绝原假设
res = sem(y,x,W);% do sem estimation
prt(res,vnames); % print the output%%结果表明lambda(which is the coefficient of the spatial dependence in the disturbances of a regression model)确实显著
%以下该检验方法lmsar想知道 在原有的SEM模型中加入spatial lag term（实质是变为类似于SAC模型）之后，spatial dependence in the residuals of 
%the sem model是否依然存在。零假设是no spatial correlation exists in the least-squares residuals of the regression model.
%检验结果表明可以拒绝零假设，就可以说 there exists spatial correlation exists in the least-squares residuals of the regression model
%after inclusion of the spatial lag term，即lambda是显著的。
res5 = lmsar(y,x,W,W);
prt(res5);
%%The first-order spatial AR model
vnames = strvcat('Duration','rho');
res = far(y,W);
prt(res,vnames);%结果表明spatial lag term的系数rho是显著的，表明不宜用不考虑the spatial lag term of dependent variable的sem模型
%%The mixed autoregressive-regressive model
%vnames = strvcat('Incidence','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','closeness','sum_Ratio09','mc-distance');
vnames = strvcat('Duration','const','Urban ratio','Density','PGDP','Income','Hospital','Bed','Doctor','College','MidSchool','PriSchool','CollegeStu','MidSchoolStu','PriSchoolStu','Closeness','sum_Ratio09','mc-distance');
res=sar(y,x,W);
prt(res,vnames);
plt(res,vnames);%结果表明加入explanatory variable之后，spatial lag of dependent variable的系数rho变得不显著了
%%SDM
res = sdm(y,x,W);
prt(res,vnames);%结果表明加入spatial lag of explanatory variable之后，spatial lag of dependent variable的系数rho依然不显著
%%SAC
W2 = slag(W,2); % standardized W2 result from slag
res1 = sac(y,x,W2,W);% general spatial model W2,W
prt(res1,vnames); % print the output
res2 = sac(y,x,W,W2);% general spatial model W,W2
prt(res2,vnames); % print the output
res3 = sac(y,x,W,W); % general spatial model W,W
prt(res3,vnames); % print the output
plt(res3);
res4 = sac(y,x,W2,W2); % general spatial model W2,W2 (selected)
prt(res4,vnames); % print the output 
