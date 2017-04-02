% Barriers
zero_1 = -3592.97102056;
E_1_l = 2.951646;
E_1_l_2 = -(-3584.167250--3581.721798);
E_1_l_0 = -3598.056400--3601.470647;
E_4_l = -3590.207247--3592.956166;
E_4_l_2 = -3581.852746--3584.077936;
E_4_l_0 = -3598.331221--3601.522984;
E_1_4 = -3589.95316882--3592.97102056;
E_3_4 = -3589.37845071-zero_1;
E_3 = -3592.97101757-zero_1;
E_3_4 = E_3_4-E_3;
E_c = -3589.86261278--3592.97102056;%3.100865 
E_1_lxp2 = -3574.65353--3577.708882;
E_1_lxm2 = -3582.279110--3585.078091;
E_1_lxm2_2 = -3573.307794--3575.558803;
E_1_lxp2_2 = -3566.993155--3569.576452;
E_1_lxp2_0 = -3582.077389--3585.555594;
E_1_lxm2_0 = -3590.988981--3594.273917;

% Boltzmann constant
kb = 8.6173324e-5;
%D = (0.0276+0.0262)/2;
%v = 109019950963228/2/pi;

% Vibrational frequency and temperature
v = 1.6e13;
T = 1300;
%T = 1800;

% Nearest neighbour distance
nnd = 3.08472680894400e-10;

% transition rates
r1 = v*exp(-E_1_l/kb/T);
r2 = v*exp(-E_4_l/kb/T);
r3 = v*exp(-E_1_4/kb/T); %4.8
r4 = v*exp(-E_3_4/kb/T);%-((E_1_l+E_4_l)/2-(4.8-3.7))/k/T);%3.7
r5 = v*exp(-E_c/kb/T);

% Silicon atom positions in unit cell
spos_Si = [0.33333333,  0.66666667,  0.93750000-1; 0.00000, 0.00000, 0.1875; 0.66666667,  0.33333333,  0.43750000; 0.000000,  0.000000,  0.68750000;  0.33333333,  0.66666667,  0.93750000; 0.00000, 0.00000, 0.1875+1];

% lattice scaling parameters and unit cell
a = 1.0104;
c = 1.0074;
cell2 = [nnd, 0*a, 0*a; -nnd/2,nnd/2*sqrt(3), 0*a; 0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220];

% Define two rotation matrix to capture all possible in-plane and out of
% plane transitions
eul1 = [pi/3,0,0];
rotm1 = eul2rotm(eul1);
eul2 = [2*pi/3,0,0];
rotm2 = eul2rotm(eul2);

% Initialize time, number of time steps (N) and number of particles (it), 
% as well as the array to store positions (llocs), which level - hexagonal
% or cubic we are on - (k) and charge state (chg)  
tt = 0;
N = 1000;
it = 1000;
k = zeros(it,1);
itt = it;
locs = zeros(1,3);
llocs = zeros(N,3,it);
% llocs(10,1:2,:) = llocs(10,1:2,:) + nnd*randn(1,2,it);
chg = -ones(1,1,it);

% Run simulation
for y = 2:N %11:N
    t = 0;
    for h = 1:it
        if k(h) == 0 || k(h) == 2
            R = 0;
            for ij = 2:7
                locs = llocs(y-1,:,h)+cell2(1,:)*rotm1^(ij-2);
                R = cat(2,R,r1*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 8:10
                locs = llocs(y-1,:,h)+(spos_Si(k(h)+1,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(ij-8);
                R = cat(2,R,r3*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 11:13
                locs = llocs(y-1,:,h)+(spos_Si(k(h)+3,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(ij-11);
                R = cat(2,R,r4*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 14:17
                R = cat(2,R,r5+R(ij-1));
            end
        elseif k(h) == 1 || k(h) == 3
            R = 0;
            for ij = 2:7
                locs = llocs(y-1,:,h)+cell2(1,:)*rotm1^(ij-2);
                R = cat(2,R,r2*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 8:10
                locs = llocs(y-1,:,h)+(spos_Si(k(h)+1,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(ij-8);
                R = cat(2,R,r4*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 11:13
                locs = llocs(y-1,:,h)+(spos_Si(k(h)+3,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(ij-11);
                R = cat(2,R,r3*exp(-(get_pot(locs,llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end]))-get_pot(llocs(y-1,:,h),llocs(y-1,:,[1:h-1,h+1:end]),chg(1,1,[1:h-1,h+1:end])))/kb/T)+R(ij-1));
            end
            for ij = 14:17
                R = cat(2,R,r5+R(ij-1));
            end
        elseif k(h) == 4
            continue
        end
        Q = R(end);
        u1 = rand;
        test = Q*u1;
        L = 0;
        E = length(R)-1;
        for j = 1:length(R)
            if L>E
                display('unsuccessful');
            end
            m = floor((E+L)/2);
            if R(m+1)<test && R(m+2)<test
                L = m+1;
            elseif R(m+1)>test && R(m+2) > test
                E = m-1;
            else
                break
            end
        end
        if m < 6
            %locs(y,:) = locs(y,:)+(locs2(y-1,:)+cell2(1,:)*rotm1^m)/it;
            llocs(y,:,h) = llocs(y-1,:,h)+cell2(1,:)*rotm1^m;
        elseif m <9
            %locs(y,:) = locs(y,:)+(locs2(y-1,:)+(spos_Si(k+1,:)-spos_Si(k+2,:))*cell2*rotm2^(m-6))/it;
            llocs(y,:,h) = llocs(y-1,:,h)+(spos_Si(k(h)+1,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(m-6);
            k(h) = mod(k(h)-1,4);
        elseif m < 12
            %locs(y,:) = locs(y,:)+(locs2(y-1,:)+(spos_Si(k+3,:)-spos_Si(k+2,:))*cell2*rotm2^(m-9))/it;
            llocs(y,:,h) = llocs(y-1,:,h)+(spos_Si(k(h)+3,:)-spos_Si(k(h)+2,:))*cell2*rotm2^(m-9);
            k(h) = mod(k(h)+1,4);
        elseif m > 11
            %locs(y,:) = locs(y,:)+(locs2(y-1,:)+(spos_Si(k+3,:)-spos_Si(k+2,:))*cell2*rotm2^(m-9))/it;
            llocs(y:end,1,h) = llocs(y-1,1,h);
            llocs(y:end,2,h) = llocs(y-1,2,h);
            llocs(y:end,3,h) = llocs(y-1,3,h);
            k(h) = 4;
            chg(1,1,h) = 2;
            itt = max([itt-1,1]);
        else
            display(['unsuccessful']);
        end
        u2 = rand;
        t = t + 1/Q*log(1/u2)/itt;
    end
    tt = tt+t;
end

% llocs = llocs*1e9;
% for y = 165%1:N-1
%     q = figure;
%     hidem(q)
%     dimenx = max(max(abs(llocs(:,1,:))));
%     dimeny = max(max(abs(llocs(:,2,:))));
%     dimenz = max(max(abs(llocs(:,3,:))));
%     dimen = max([dimenx,dimeny,dimenz]);
%     for l = 1:it
%         if llocs(y,:,l) == llocs(y+1,:,l) & y >1
%             scatter3(llocs(y,1,l),llocs(y,2,l),llocs(y,3,l),'r.')
%         else
%             scatter3(llocs(y,1,l),llocs(y,2,l),llocs(y,3,l),'b.')
%         end
%         hold on
%     end
%     %quiver3(locs(1,1),locs(1,2),locs(1,3),locs(y+1,1),locs(y+1,2),locs(y+1,3))
% %     xlim([-sqrt(N)*cell2(1,1)/sqrt(it),sqrt(N)*cell2(1,1)/sqrt(it)])
% %     ylim([-sqrt(N)*cell2(1,1)/sqrt(it),sqrt(N)*cell2(1,1)/sqrt(it)])
% %     zlim([-sqrt(N)*cell2(1,1)/sqrt(it),sqrt(N)*cell2(1,1)/sqrt(it)])
%     xlim(1*[-dimen,dimen])
%     ylim(1*[-dimen,dimen])
%     zlim(1*[-dimen,dimen])
%     xlabel('\Deltax (nm)','Fontsize',16)
%     ylabel('\Deltay (nm)','Fontsize',16)
%     zlabel('\Deltaz (nm)','Fontsize',16)
%     title(['after ',num2str(y/N*tt,'%.2f'),' seconds'],'Fontsize',18)
%     set(gcf,'Color','w')
%     axis square
%     saveas(gcf,['Video/im',num2str(y,'%03.0f'),'.jpg'])
%     delete(q)
% end    


%E_Siv0 = 4.1;
%E_Siv2 = 3.6;
%E_Siv1 = (E_Siv0+E_Siv2)/2;



    