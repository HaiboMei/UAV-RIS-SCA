%configure GTs' distribution.
%most of the code is from the LTE simulator;
function [W_k,K] =GT_distributions(distance)  
    K=0;
    if (1)  %GTs in Cellular based distribuion.  
        inter_bts_distance = distance/2;
        GTs_per_cell=6;
            
        % Number of BTS rings in the network map
        n_rings = 1;
        [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
            (-n_rings:n_rings)*sin(pi/3)); %regular grid
        if mod(n_rings,2) == 0
            tmp_shift_idx = 2:2:2*n_rings+1; %shift all even rows
        else
            tmp_shift_idx = 1:2:2*n_rings+1; %shift all odd rows
        end
        tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5; %shift

        rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; %rotation operator
        for i_ = 1:7
            %border of the network
            tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).';
        end
        
        tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));
        tmp_x = tmp_gridx(tmp_valid_positions);
        tmp_y = tmp_gridy(tmp_valid_positions);
        
        eNodeB_positions = [distance+tmp_x*inter_bts_distance distance+tmp_y*inter_bts_distance];  
                
        for e=1:size(eNodeB_positions,1)
            eNodeB_positions(e,1) = eNodeB_positions(e,1)+randi([-10,10]);
            eNodeB_positions(e,2) = eNodeB_positions(e,2)+randi([-10,10]);
        end
        
        K=(size(eNodeB_positions,1)-4)*GTs_per_cell;                   %Number of GTs
        W_k=zeros(2,K); 
        radius=10;
        for e=1:(size(eNodeB_positions,1)-4)
            for k=1:GTs_per_cell
                W_k(1,(e-1)*GTs_per_cell+k) = eNodeB_positions(e,1)+randi([-radius,radius]);
                W_k(2,(e-1)*GTs_per_cell+k) = eNodeB_positions(e,2)+randi([-radius,radius]);
            end
        end
        %plot the GT distributions for anlysis purpose.
    end
    if (0)
        K=6;                   %Number of GTs
        W_k=zeros(2,K); 
        for k=1:K
            W_k(1,k) = randi([0,distance]);
            W_k(2,k) = randi([0,distance]);
        end
    end
end