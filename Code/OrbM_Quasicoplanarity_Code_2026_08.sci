//Scilab code for Computational Analysis of Quasi-planarity of Orbital 
//Trajectories of Marbles on a 3D Printed Inverse-r Surface.
//----------------------------------CODE----------------------------------
clear
clc

codeDir = get_absolute_file_path("OrbM_Quasicoplanarity_Code_2026_08.sci");
dataFile = codeDir + "../Tracked_Orbit_Reference_Set/Closed_Orbit_30_FPS.xls";
data = readxls(dataFile);

s0 = data(1)
T_1 = s0.value()
T_2 = s0.text()
Track_Data = T_1(3:$, 2:3)
xT = Track_Data(:,1)
yT = Track_Data(:,2)
N = length(xT)

c1 = 1.0648
c2 = -12.071687
c3 = 6.42085

for i = 1:1:N
    rT(i) = sqrt(xT(i)^2 + yT(i)^2)
    zT(i) = c1 + c2/rT(i)
end

// --------------------------------------------------------
// ------ Angles Assessment, and Quasicoplanarity ---------
for j = 1:1:N-2
    rf1 = [xT(j) yT(j) zT(j)]
    rf2 = [xT(j+1) yT(j+1) zT(j+1)]
    rf3 = [xT(j+2) yT(j+2) zT(j+2)]
    
    rf23 = rf3 - rf2
    rf12 = rf2 - rf1
    nf = cross(rf12, rf23)
    unf = nf/norm(nf)
    CrossVectors(:,j) = unf
    N_Matrix(j,:)=unf
end

//Angles Assessment - Method 1: Angles between consecutive normal vectors
for i = 1:N-3
    tv = CrossVectors(:,i)'*CrossVectors(:,i+1);
    Angles_Method1(i) = acos(tv)*180/%pi;
end
Index_N3 = 1:N-3
Mean_Angles_Method1 = mean(Angles_Method1)
Sig_Method1 = stdev(Angles_Method1)//sample stdev
dof = length(Angles_Method1)-1; //deg of freedom = N3-1
Std_Err_Theta_Method1 = Sig_Method1/sqrt(length(Angles_Method1)); //standard error in mean, standard deviation divided by square root of N
disp('-----------------------')
disp('Students t-test for Angles Assessment:')
disp('Method 1 (Angle b/w Consecutive Normals)') 
disp('(Sample) Mean Angle = ', Mean_Angles_Method1)
disp('(Sample) Standard Deviation = ', Sig_Method1)
disp('Standard Error in mean angle = ', Std_Err_Theta_Method1)
cutoff_angle = 10;//theta_0 in t test
disp('Degrees of Freedom', dof)
t_statistic = (Mean_Angles_Method1 -cutoff_angle)/Std_Err_Theta_Method1;
disp('t-statistic = ', t_statistic)

[tP,tQ] = cdft("PQ",t_statistic,dof)
p_value = tP; //lower tailed comparison
disp('p-value = ',p_value)
Conf_Lev = 0.999;
Conf_Lev_Per = 100*Conf_Lev;
alpha_value = 1-Conf_Lev;

if (p_value<alpha_value), disp('We reject the null hypothesis.'),
    disp ('These normal vectors can be regarded as parallel,') 
    disp('and points as coplanar, with confidence level (%) = ', Conf_Lev_Per) 
else  disp('We fail to reject the null hypothesis.') 
    disp('Insuffient statistical evidence to suggest')
    disp('that the normals as parallel, and so, points')
    disp('may not be regarded as approximately coplanar')
end

subplot(1,2,1)
plot(Index_N3,Angles_Method1,'o')
plot(Index_N3,Mean_Angles_Method1,'red')
subplot(1,2,2)
histplot(10,Angles_Method1, normalization=%f)

//Angles Assessment - Method 2: Angles between every pair of normal vectors
cc = 0
for kk = 1:N-3
    for kkk = kk+1:N-2
       cc = cc+1
       tv3 = acos(abs(CrossVectors(:,kk)'*CrossVectors(:,kkk)))*180/%pi;
       Angles_Every_Pair(cc) = tv3;
       Angles_Pair_Matrix_Form(kk,kkk)=tv3;
    end
end
Index_3 = 1:cc
Mean_Angles_Method2 = mean(Angles_Every_Pair);
figure
subplot(1,2,1)
plot(Index_3,Angles_Every_Pair,'o')
ff = gcf()
ff.background = -2;
plot(Index_3,Mean_Angles_Method2,'red')
subplot(1,2,2)
histplot(10,Angles_Every_Pair, normalization=%f)
Mean_Angles_Method2 = mean(Angles_Every_Pair);
Sig_Method2 = stdev(Angles_Every_Pair)//sample stdev
dof = length(Angles_Every_Pair)-1; //deg of freedom = N3-1
Std_Err_Theta_Method2 = Sig_Method2/sqrt(length(Angles_Every_Pair)); //standard error in mean, standard deviation divided by square root of N
disp('-----------------------')
disp('Students t-test for Angles Assessment:')
disp('Method 2 (Angle b/w every pair of normals)') 
disp('Sample size, n = ',cc)
disp('(Sample) Mean Angle = ', Mean_Angles_Method2)
disp('(Sample) Standard Deviation = ', Sig_Method2)
disp('Standard Error in mean angle = ', Std_Err_Theta_Method2)
cutoff_angle = 10;//theta_0 in t test
disp('Degrees of Freedom', dof)
t_statistic = (Mean_Angles_Method2 -cutoff_angle)/Std_Err_Theta_Method2;
disp('t-statistic = ', t_statistic)

[tP,tQ] = cdft("PQ",t_statistic,dof)
p_value = tP; //lower tailed comparison
disp('p-value = ',p_value)
Conf_Lev = 0.999;
Conf_Lev_Per = 100*Conf_Lev;
alpha_value = 1-Conf_Lev;

if (p_value<alpha_value), disp('We reject the null hypothesis.'),
    disp ('These normal vectors can be regarded as parallel,') 
    disp('and points as coplanar, with confidence level (%) = ', Conf_Lev_Per) 
else  disp('We fail to reject the null hypothesis.') 
    disp('Insuffient statistical evidence to suggest')
    disp('that the normals as parallel, and so, points')
    disp('may not be regarded as approximately coplanar')
end
disp('-----------------------')

disp("Z-test for Angles Assessment")
disp('Method 2 (Angle b/w every pair of normals)') 
disp('Sample size, n = ',cc)
disp('(Sample) Mean Angle = ', Mean_Angles_Method2)
disp('(Sample) Standard Deviation = ', Sig_Method2)
disp('Standard Error in mean angle = ', Std_Err_Theta_Method2)
sig_theta = stdev(Angles_Every_Pair);
disp('Standard Deviation of the Angle between Normals',sig_theta)
disp('So, as per this calculation, s = sig_theta/sqrt{n} =', sig_theta/cc^0.5)
variance_theta = variance(Angles_Every_Pair);
s_stddev = (variance_theta./length(Angles_Every_Pair)).^0.5;
disp('Alternative calculation gives s = sigma/sqrt{N} = ', s_stddev)

z_statistic = (Mean_Angles_Method2 -cutoff_angle)/s_stddev;
disp('The Z-statistic is calculated to be:', z_statistic)

[zP,zQ]=cdfnor("PQ",z_statistic,0,1);//Against N(0,1) normal distribution
p_value_Z = zP; //lower tailed comparison
disp('p-value (X 10^-7) = ',p_value_Z*1D7)
Conf_Lev = 0.999;
Conf_Lev_Per = 100*Conf_Lev;
alpha_value = 1-Conf_Lev;

if (p_value_Z<alpha_value), disp('We reject the null hypothesis.'),
    disp ('These normal vectors can be regarded as parallel,') 
    disp('and points as coplanar, with confidence level (%) = ', Conf_Lev_Per) 
else  disp('We fail to reject the null hypothesis.') 
    disp('Insuffient statistical evidence to suggest')
    disp('that the normals as parallel, and so, points')
    disp('may not be regarded as approximately coplanar')
end
disp('-----------------------')
//-----------------------------------------------------------------------------

// ---------------------OTHER FIGURES------------------------------ //
//xdel(winsid())
figure
ff = gcf()
ff.background = -2;
//An animation of the mapped/recovered 3D coordinates of marble motion.
second = 3
comet3d(xT, yT, zT)
fig = gca()
fig.box = "on"
fig.cube_scaling = "on"
fig.rotation_angles = [35,45]
fig.children(1).mark_size = 15
fig.children(1).mark_background = -2
fig.children(2).mark_mode = "on"
fig.children(2).mark_style = 0
fig.children(2).mark_size_unit = "point"
fig.children(2).mark_size = 4
fig.children(2).mark_foreground = 5
fig.children(2).mark_background = 5
fig.data_bounds=[-8,-9.5,-2;6,6,0.01];
fig.sub_ticks = [0,0,0]
sleep(second,"s")

figure
ff = gcf()
ff.background = -2;
//Marble points atop the surface
scatter3d(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")
fig = gca()
fig.box = "on"
fig.cube_scaling = "on"
fig.rotation_angles = [35,45]
fig.children.mark_style = 9
fig.children.mark_size_unit = "point"
fig.children.mark_size = 6
fig.data_bounds=[-8,-9.5,-2;6,6,0.01];
fig.sub_ticks = [0,0,0]
xgrid(8,0,10)

//Adding in the surface too
xx = linspace(-10,10,51);
yy = linspace(-10,10,51);
for i = 1:1:length(xx)
    for j = 1:1:length(yy)
    rr(i,j) = sqrt(xx(i)^2 + yy(j)^2)
    zz(i,j) = c1 + c2/rr(i,j)
end
end

surf(xx,yy,zz)
xgrid(8,0,10)
bbb = gca()
bbb.box = "back_half"
bbb.tight_limits = ["on","on","on"]
bbb.data_bounds(5) = -10 //stop z at -10 only
bbb.grid = [0,0,0]
aaar = bbb.children(1)
aaar.foreground = 16
aaar.color_mode = 7
aaar.color_flag = 0
sleep(second,"s")

figure
ff = gcf()
ff.background = -2;
//A static sketch containing the 3D points at different times, in space (3D figure).
scatter3d(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")
fig = gca()
fig.box = "on"
fig.cube_scaling = "on"
fig.rotation_angles = [35,45]
fig.children.mark_style = 9
fig.children.mark_size_unit = "point"
fig.children.mark_size = 4
fig.data_bounds=[-8,-9.5,-2;6,6,0.01];
fig.sub_ticks = [0,0,0]
xgrid(8,0,10)
scatter3d(xT, yT, zT)
fig2 = gca()
aab = fig2.children(1)
aab.mark_style = 9;
aab.mark_size = 12;
aab.mark_foreground = 2;
aab.mark_foreground = 2;

sleep(second,"s")

figure
ff = gcf()
ff.background = -2;

scatter3d(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")

scatter3d(xT,yT,zT) 
//One for the red marble center dots, and second for the enclosing blue circle - finite marble size


kk = 0:0.1:1;

xN = zeros(N-2,11);
yN = zeros(N-2,11);
zN = zeros(N-2,11);

for i = 1:N-2
    tv = N_Matrix(i,:); //temporary variable
    xN(i,:) = xT(i+1)-tv(1)*kk;
    yN(i,:) = yT(i+1)-tv(2)*kk;
    zN(i,:) = zT(i+1)-tv(3)*kk;
end

//Normal vector lines
scatter3d(xN,yN,zN)

//Endpoints of the normal vectors - For the upward facing arrows
for i = 1:N-2
    xEnd(i) = xN(i,$);
    yEnd(i) = yN(i,$);
    zEnd(i) = zN(i,$);
end

scatter3d(xEnd,yEnd,zEnd)

aaa = gca()

aaaa = aaa.children(1)
aaaa.mark_style = 6;
aaaa.mark_size = 9;
aaaa.mark_foreground = 2;
aaaa.mark_foreground = 2;

aaab = aaa.children(2)
aaab.mark_style = 0;
aaab.mark_size = 3;
aaab.mark_foreground = 2;
aaab.mark_foreground = 2;

aaac = aaa.children(3)
aaac.mark_size = 10;

aaad = aaa.children(4)
aaad.foreground = 5;
aaad.line_mode = 'ON';
aaad.line_style = 2;
aaad.mark_size = 3;

