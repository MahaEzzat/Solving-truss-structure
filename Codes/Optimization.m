clear
clc
syms R1 R2 R17 
h_range = 0.0001:0.0001:0.0068;
A_range = 0.00001:0.00001:0.001;

for count_A=1:length(A_range) 
for count_h=1:length(h_range)
%%Assumed
h = h_range(count_h) %Beam height in meter
A3 = A_range(count_A)     %element (8 9 12 13) area in m^2
A1 = A3;           %element (6 7 14 45 16 19) area in m^2
A2 = A3;           %element (17 18) area in m^2
A4 = 0.00000001*A3; %element (10 11) area Zero force members in m^2 (almost zero)


Ab = 3.5*h; %meter^2
Ib=(3.5*h^3)/12; %meter^4
Beam_density = 2300; %kg/m^3
Truss_density = 7800; %kg/m^3


% Young’s modulus for the element material 
Eb=41*10^9; %Beam  E
Et=207*10^9;%Truss E
le=2; %meter
w=-0.9*9.81*1000; %N/m

%% L global
%The angle corresponding to each element 
theta=[0 0 0 0 0 120 60 120 60 120 60 120 60 120 60 0 0 0 0 ];
lb =[0 w*le/2 (w*le^2)/12 ]';
lb2=[0 w*le/2 -(w*le^2)/12 ]';
R=[R1;R2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;R17;0];
Lb=[lb;lb2+lb;lb2+lb;lb2+lb;lb2+lb;lb2]+R;
Lt=[0 0 0 0 0 0 0 0 0 0]';
L=[Lb;Lt];

%% local K's
for i=1:19
    if i<6
     m=Eb*Ib/le^3;
     n=Ab*Eb/le;
    K{i}=[n 0 0 -n 0 0 ;...
        0 m*12 m*6*le 0 -12*m 6*le*m;...
        0  6*m*le 4*m*le^2 0 -6*m*le 2*m*le^2;...
        -n 0 0 n 0 0 ; ...
        0  -12*m -6*m*le 0  12*m -6*m*le;...
        0 6*m*le 2*m*le^2 0  -6*m*le 4*m*le^2];% stiffness matrix of the beam
    else
        %changing area
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
    K{i}=(A*Et/le)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;... % stiffness matrix of the truss
   -C*C -C*S C*C C*S; -C*S -S*S C*S S*S];
    end
end
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

%% row & column vectors
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

%% K global
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
u=vpa(u);%deltas


%% Internal forces in global coordinates & weight
 for ll=1:19 %for each element
    C=cosd(theta(ll));
    S=sind(theta(ll)); 
    T=[C S 0 0;...
        -S C 0 0;...
        0 0 C S;...
        0 0 -S C];
    
      %changing area
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
    if sigma_beam(ll)> 1.3*10^8
        break
    end
    Weight(ll) = Beam_density*Ab*le;
    
     elseif ll>=6 && ll<16
    uu2(:,ll-5)=u(c2(ll-5,:));
    f2(:,ll-5)=K{ll}*uu2(:,ll-5);
    f2_axial(:,ll-5)=T*f2(:,ll-5);
    sigma_truss(ll-5) = abs(f2_axial(1,ll-5))/A;
    
    if sigma_truss(ll-5)> 2.2*10^8
        break
    end
     
    Weight(ll) = Truss_density*A*le;
    
     elseif ll>=16
    uu2(:,ll-5)=u(c3(ll-15,:)) ; 
    f2(:,ll-5)=K{ll}*uu2(:,ll-5);% internal forces for trusses
    f2_axial(:,ll-5)=T*f2(:,ll-5);
    sigma_truss(ll-5) = abs(f2_axial(1,ll-5))/A;
    
     if sigma_truss(ll-5)> 2.2*10^8
        break
     end
     
    Weight(ll) = Truss_density*A*le;
    
    end
 end
 
  if ll==19 
 f_beam = double(f); %beam internal forces
 f_axial_truss=double(f2_axial);%axial forces in the truss members
 sigma_beam_total(count_h,:,count_A) = double(sigma_beam); %Beam Stress in N/m2
 sigma_truss_total(count_h,:,count_A) = double(sigma_truss); %Truss Stress in N/m2
 Total_Weight(count_h,count_A) = sum(Weight) %kg
  end
  
end
end

%% Results
Min_weight = min(Total_Weight(Total_Weight>0))
[row,col] = find(Total_Weight==Min_weight);
A = A_range(col)
h = h_range(row)