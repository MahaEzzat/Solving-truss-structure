clear
clc
syms R1 R2 R17  %Reactions of pinned & roller supports
syms e   
 
%% Assumptions
A3 = 2.853400e-04;  %Element (17 18) area in m^2
A1 = 0.6838*A3;     %Element (6 7 14 45 16 19) area in m^2
A2 = 0.31621*A3;    %Element (8 9 12 13) area in m^2
A4 = 0.00000001*A3; %Element (10 11) area Zero force members in m^2 (almost zero)
h = 0.0062532;      %Beam height in meter

%% Givens
Ab = 3.5*h;           %meter^2
Ib=(3.5*h^3)/12;      %meter^4
Beam_density = 2300;  %kg/m^3
Truss_density = 7800; %kg/m^3

Eb=41*10^9;       %Beam Young’s modulus
Et=207*10^9;      %Truss Young’s modulus
le=2;             %Element length in meter
w=-0.9*9.81*1000; %Distributed load in N/m

%% L global
theta=[0 0 0 0 0 120 60 120 60 120 60 120 60 120 60 0 0 0 0 ]; %The angle corresponding to each element 
lb =[0 w*le/2 w*le^2/12 ]';                  %Distributed load contribution
lb2=[0 w*le/2 -w*le^2/12 ]';                 %Distributed load contribution
R=[R1;R2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;R17;0]; %Reactions contribution
Lb=[lb;lb2+lb;lb2+lb;lb2+lb;lb2+lb;lb2]+R;   %Beam Load Matrix
Lt=[0 0 0 0 0 0 0 0 0 0]';                   %Truss Load Matrix
L=[Lb;Lt];                                   %Global Load Matrix

%% local stiffness matrices
for i=1:19
    if i<6
     m=Eb*Ib/le^3;
     n=Ab*Eb/le;
    K{i}=[n 0 0 -n 0 0 ;...
        0 m*12 m*6*le 0 -12*m 6*le*m;...
        0  6*m*le 4*m*le^2 0 -6*m*le 2*m*le^2;...
        -n 0 0 n 0 0 ; ...
        0  -12*m -6*m*le 0  12*m -6*m*le;...
        0 6*m*le 2*m*le^2 0  -6*m*le 4*m*le^2];   %stiffness matrix of the beam
    else
        
        %changing truss area
      if (i == 6 || i == 7 || i == 14 || i == 15 || i == 16 || i == 19)
          A = A1;
      end
       if (i == 8 || i == 9 || i == 12 || i == 13)
          A = A2;
       end
       if (i == 17 || i == 18)
          A = A3;
       end
       if (i == 10 || i == 11)
          A = A4;
       end
      
    C=cosd(theta(i));
    S=sind(theta(i));
    K{i}=(A*Et/le)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;... %stiffness matrix of the truss
   -C*C -C*S C*C C*S; -C*S -S*S C*S S*S];
    end
end

%stiffness matrix in vector form
iii=1;
for ii=1:19
for jj=1:length(K{ii}) % raw loop
for iij=1:length(K{ii}) %col. loop
    
    p(iii,1)=K{ii}(jj,iij);
iii=iii+1;
end
end
end
Kvec=p; 

%% Row & column vectors
c1=[1 2 3 4 5 6;4 5 6 7 8 9;7 8 9 10 11 12;10 11 12 13 14 15;13 14 15 16 17 18];
i=0;
%Beam elements 1<5
for nb=1:5;
for ndof=1:6;
    r1(1+6*i)=c1(nb,ndof);
    r1(2+6*i)=c1(nb,ndof);
    r1(3+6*i)=c1(nb,ndof);
    r1(4+6*i)=c1(nb,ndof);
    r1(5+6*i)=c1(nb,ndof);
    r1(6+6*i)=c1(nb,ndof);
    
    col1(1+6*i)=c1(nb,1);
    col1(2+6*i)=c1(nb,2);
    col1(3+6*i)=c1(nb,3);
    col1(4+6*i)=c1(nb,4);
    col1(5+6*i)=c1(nb,5);
    col1(6+6*i)=c1(nb,6);
    i=i+1;
end
end

%Truss elements from 6-->15
c2=[19 20 1 2; 19 20 4 5; 21 22 4 5;21 22 7 8; 23 24 7 8; 23 24 10 11; 25 26 10 11; 25 26 13 14; 27 28 13 14; 27 28 16 17];
ii=0;
for nt=1:10
for ntdof=1:4
    r2(1+4*ii)=c2(nt,ntdof);
    r2(2+4*ii)=c2(nt,ntdof);
    r2(3+4*ii)=c2(nt,ntdof);
    r2(4+4*ii)=c2(nt,ntdof);
    
    col2(1+4*ii)=c2(nt,1);
    col2(2+4*ii)=c2(nt,2);
    col2(3+4*ii)=c2(nt,3);
    col2(4+4*ii)=c2(nt,4);
    ii=ii+1;
end
end

%Truss elements from 16-->19
c3=[19 20 21 22; 21 22 23 24; 23 24 25 26; 25 26 27 28];
iii=0;
for nt2=1:4
for nt2dof=1:4
    r3(1+4*iii)=c3(nt2,nt2dof);
    r3(2+4*iii)=c3(nt2,nt2dof);
    r3(3+4*iii)=c3(nt2,nt2dof);
    r3(4+4*iii)=c3(nt2,nt2dof);
    
    col3(1+4*iii)=c3(nt2,1);
    col3(2+4*iii)=c3(nt2,2);
    col3(3+4*iii)=c3(nt2,3);
    col3(4+4*iii)=c3(nt2,4);
    iii=iii+1;
end
end
r=[r1 r2 r3]';
col=[col1 col2 col3]';

%% Global stiffness matrices
s=sparse(r,col,Kvec);
k_global= full(s);

k_global(1,:)=[]; k_global(:,1)=[]; L(1)=[];
k_global(1,:)=[]; k_global(:,1)=[]; L(1)=[];
k_global(15,:)=[]; k_global(:,15)=[]; L(15)=[];

%% Deltas
u=k_global\L;
u = [0; 0; u];
u1=u(1:16);
u2=u(17:27);
u=[u1;0;u2];
u=double(u);%global deflections solution


%% Internal forces in global coordinates, stresses & weights
 for ll=1:19 %for each element
    C=cosd(theta(ll));
    S=sind(theta(ll)); 
    T=[C S 0 0;...
        -S C 0 0;...
        0 0 C S;...
        0 0 -S C];
    
      %changing truss area
      if (ll == 6 || ll == 7 || ll == 14 || ll == 15 || ll == 16 || ll == 19)
          A = A1;
      end
       if (ll == 8 || ll == 9 || ll == 12 || ll == 13)
          A = A2;
       end
       if (ll == 17 || ll == 18)
          A = A3;
       end
       if (ll == 10 || ll == 11)
          A = A4;
       end 

     if ll<6
    uu(:,ll)=u(c1(ll,:));
    f(:,ll)=K{ll}*uu(:,ll);% internal forces for beams w.r.t global coordinates
    sigma_beam(ll) = abs(f(1,ll))/Ab + (max([abs(f(3,ll)) abs(f(6,ll))])*h/2)/Ib;
    Weight(ll) = Beam_density*Ab*le;
    
%% Beam Deflection function
Vi = uu(2,ll); Vj = uu(5,ll); thetai = uu(3,ll); thetaj = uu(6,ll); ui = uu(1,ll); uj = uu(4,ll);
 H1 = 0.5 - 3*e/4 + (e^3)/4;
 H3 = 0.5 + 3*e/4 - (e^3)/4;
 H2 = (1/8)*( 1 - e - e^2 + e^3);
 H4 = (1/8)*(-1 - e + e^2 + e^3);
 N1 = (1-e)/2;
 N2 = (1+e)/2;
 Y_def_function(ll) = vpa(Vi*H1 + Vj*H3 + thetai*H2 + thetaj*H4);
 X_def_function(ll) = vpa(ui*N1 + uj*N2);
 Theta_function(ll) = (2/le)*diff(Y_def_function(ll),e);
 M_function(ll) = (Eb*Ib*2/le)*diff(Theta_function(ll),e);  
 V_function(ll) = (2/le)*diff(M_function(ll),e);
 
     elseif ll>=6 && ll<16
    uu2(:,ll-5)=u(c2(ll-5,:));
    f2(:,ll-5)=K{ll}*uu2(:,ll-5);
    f2_axial(:,ll-5)=T*f2(:,ll-5);
    sigma_truss(ll-5) = abs(f2_axial(1,ll-5))/A;
    Weight(ll) = Truss_density*A*le;
    
     elseif ll>=16
    uu2(:,ll-5)=u(c3(ll-15,:)) ; 
    f2(:,ll-5)=K{ll}*uu2(:,ll-5);% internal forces for trusses
    f2_axial(:,ll-5)=T*f2(:,ll-5);
    sigma_truss(ll-5) = abs(f2_axial(1,ll-5))/A;
    Weight(ll) = Truss_density*A*le;
     end
 end
 
 %% Ploting
 for i = 1:11 %nodes position after deflections
     if i<7
 X(i) = (i-1)*le;                                %Intail position
 Y(i) = 0;                                       %Intail position       
 X_deflection(i) = u(1+3*(i-1))+(i-1)*le;        %Deflection + Intail position
 Y_deflection(i) = u(2+3*(i-1));                 %Deflection + Intail position
 Theta_deflection(i) = u(3+3*(i-1));
     else
 X(i) = (9-(i-7)*le);                            %Intail position
 Y(i) = -(le/sind(60));                          %Intail position      
 X_deflection(i) = u(27-2*(i-7))+(9-(i-7)*le);   %Deflection + Intail position
 Y_deflection(i) = u(28-2*(i-7))-(le/sind(60));  %Deflection + Intail position
     end
 end
 
 for i = 1:6
   X(2*(i-1)+12) = X(7-i);                       %Intail position
   Y(2*(i-1)+12) = Y(7-i);                       %Intail position
   X_deflection(2*(i-1)+12) = X_deflection(7-i); %Deflection + Intail position
   Y_deflection(2*(i-1)+12) = Y_deflection(7-i); %Deflection + Intail position
   if i~=6
   X(2*(i-1)+13) = X(6+i);                       %Intail position
   Y(2*(i-1)+13) = Y(6+i);                       %Intail position
   X_deflection(2*(i-1)+13) = X_deflection(6+i); %Deflection + Intail position
   Y_deflection(2*(i-1)+13) = Y_deflection(6+i); %Deflection + Intail position
   %Beam Deflections Plotting
   plot( (i-1)*le:0.1:i*le ,subs(Y_def_function(i),e,-1:0.1:1),'color',[0.2 0 0.5]);
   title('Beam Y_ _d_e_f_l_e_c_t_i_o_n Vs. X');
   ylabel('Y_ _d_e_f_l_e_c_t_i_o_n');
   xlabel('X_a_x_i_s');
   hold on
   end
 end
 
 %Beam Theta Defleaction Plotting
 figure
 for i = 1:5
   plot( (i-1)*le:0.1:i*le ,subs(Theta_function(i),e,-1:0.1:1),'color',[0.8 0.2 0.1]);
   title('Beam Theta_ _d_e_f_l_e_c_t_i_o_n Vs. X');
   ylabel('Theta_ _d_e_f_l_e_c_t_i_o_n');
   xlabel('X_a_x_i_s');
   hold on    
 end
 
 %Beam Bending Moment Plotting
 figure
 for i = 1:5
   plot( (i-1)*le:0.1:i*le ,subs(M_function(i),e,-1:0.1:1),'color',[0.8 0.8 0.1]);
   title('Beam Bending Moment Vs. X');
   ylabel('Bending Moment_(_N_._m_)');
   xlabel('X_a_x_i_s');
   hold on    
 end
 
  %Beam Shear force Plotting
 figure
 for i = 1:5
   plot( (i-1)*le:0.1:i*le ,subs(V_function(i),e,-1:0.1:1),'color',[1 0.1 0.1]);
   title('Beam Shear Force Vs. X');
   ylabel('Shear Force_(_N_)');
   xlabel('X_a_x_i_s');
   hold on    
 end
 
 
 
 
 %Truss before Deflections
 figure
 plot(X(1:11),Y(1:11),'color',[0 0 1]);
 hold on
 plot(X(12:22),Y(12:22),'color',[0 0 1]);
 title('Truss before Deflections');
 ylabel('Y-axis');
 xlabel('X-axis');
 hold on
 %Truss after Deflections
 plot(X_deflection(1:11),Y_deflection(1:11),'color',[1 0 0]);
 hold on
 plot(X_deflection(12:22),Y_deflection(12:22),'color',[1 0 0]);
 title('Truss after Deflections');
 ylabel('Y-axis');
 xlabel('X-axis');
 ylim([-2.4 0.1]);
 
 for i=1:11
     if i<7
    text(X_deflection(i),Y_deflection(i), ['Node '  num2str(i)],'FontSize',8);
     else
     text(X_deflection(i),Y_deflection(i), ['Node '  num2str(18-i)],'FontSize',8);   
     end
 end
 
 %% Results
 f_beam = double(f) %beam internal forces
 f_axial_truss=double(f2_axial)%axial forces in the truss members
 sigma_beam = double(sigma_beam) %Beam Stress in N/m2
 sigma_truss = double(sigma_truss) %Truss Stress in N/m2
 Total_Weight = sum(Weight) %kg
 

 