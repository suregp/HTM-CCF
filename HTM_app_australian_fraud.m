%% NOTE******:%This program is strictly for research purposes and should be used 
%with care. The Authors will not provide any warranty for resulting
%damages from use of this software


%Author: E.N. Osegi
%Affiliation: National Open University of Nigeria(NOUN)
%Title: HTM SPATIAL-POOLER ALGORITHM FOR CREDIT CARD FRAUD (CCF) DETECTION
%Version: v1
%Date: 15th June, 2021


%%
clc;
clear all;

per_cent = 70;%

%% INPUT: data:

data_name = 'australian_data_full_690'; %Note 690 samples
%% HTM-SP PARAMETERS:
iters = 5;% Values: 5,10,30,100,1000 (depend on computational capacity)
initPerm = 0.21; %Initial Permanence
connectedPerm_ = 0.51; % Permanence threshold (Key - improper setting may cause error)
synPermActiveInc = 0.1; %Synapse Permanence Incrementt
synPermInactiveDec = 0.02; %Synapse Permanence Decrementt
n = 256; % Consider 256, 512, 1024, 2048 Typical
w = 2; % sdr vector number of ON bits (Key - improper setting may cause error)
%Q_th = w; %stimulus threshold
Q_th = w; % Threshold for inhibition (winner columns) is
          %made == w for exact matches
          % For Inexact matches(more practical), Qth < w e.g Qth = w-1 if
          % w = 2 , or Qth = w-2 if w =3 or in any suitable fashion;
boost_factor = 100; % Start-up value for the boost
%Note: s = w/n;
col_size = n; %Should be >> 250; Typical vals: 250, 500, 1000, 2048: see pg 49 of 69 of BAMI
stimulusThreshold  = Q_th; % typical 0 to 5 ; see pg 49 of 69 of BAMI
desired_localActivity = roundn(0.02*col_size,1); %see pg 49 of 69 of BAMI book
                                                  %desired_localActivity ==
                                                  %numActiveColumnsPerInhAr
                                                  %ea 

a = 0.0; b = 1.0; %Random width
activeColumns = zeros(1,col_size);
activeColumns_z = zeros(iters,col_size);   
boost = boost_factor.*(ones(1,col_size));
learn_state = 1; %Set to 0 to turn-off learn state                                                  



%% Data Handling:
data_name_n =  [data_name '.' 'txt'];

B =  textread(data_name_n, '%s', 'delimiter', '\n', ...
                'whitespace', '');

input_data_n = B;%cellstr(num2str(roundn((xlsread('data_in.xls')),-1)));
%input_data_n = {1,2,3};% 

len_input_data = length(input_data_n);

samples_n = roundn((per_cent/100)*len_input_data,0);
%input_data = input_data_n(t_o);
%input_data_k = input_data(t_o);

%% INPUT:
%%
l_c_range = 1;
[roBo,coBo] = size(B);
l_c_rangeo = roBo; %

%APstoredevo = zeros(roBo,l_c_rangeo);
for to_no = 1:roBo

 kodevo =  double(cell2mat(B(to_no,:)));
 lo_apno(to_no) = length(kodevo);
 APstoredevo(to_no,1:lo_apno(to_no)) =  kodevo;
 %string_label(to_no,:) = (cell2mat(input_data_n(to_no,:)));
 %class_label(to_no,1) = string_label(to_no,lo_apno(to_no));
 class_label(to_no,1) = ((APstoredevo(to_no,lo_apno(to_no)))-1)==49;

end


%% Encoding:
%% Transformation to Binary Chains:    
%c_mat = logical(reshape_nd(st_bin_val,1,ro_bin*co_bin))';
%Normalization:
APstoredevo_prob = APstoredevo/127; %127 --> the maxima of the ASCII Chart
uoo = APstoredevo./max(max(APstoredevo));

%for j = 1:roBo
j = roBo;
    for i = 1:roBo
        A_j(i,:) = j;
        A_i(i,:) = i;
        k1(i,:) = uoo(j,:);
        k2(i,:) = uoo(i,:);
        Agn(i,:) = (uoo(j,:)==uoo(i,:));%% j == nth observation
                                        %% also you can set j == 1st
                                        %% observation during the encoding
                                        %% run
    end
%end
Agn = [Agn class_label];% Super-imposition of class labels
%%

%% SPATIAL-POOLING BEGINS HERE:
for t_o = 1:l_c_rangeo 


    if(t_o <= samples_n)
        %c_mato = APstoredevo_bin_transpose(:,t_o);
        c_mato = Agn(t_o,:)';
        %[ro_col, co_col] = size(c_mato);
        [ro_col, co_col] = size(c_mato);
        
        %% Control signal for learning:
        if(learn_state == 1)
            iters = iters+1;
        else
            iters = 1;
        end
       
                
        %% Training:
        for it = 1:iters
            %% Initializing Permanences:
            if(it > 1)
                connectedPerm = connectedPerm_;
            else
                connectedPerm = initPerm;
            end
            
        %% Overlap:                                                 


                    for c = 1:col_size
                        %% Initialize Permanences:    
                        Overlap_z = 0; % Overlap-Generator Initialization 
                        permanence(:,c) = (a + (b-a) * rand(ro_col,1));
                        synapse.perms(:,c) = (a + (b-a) * rand(ro_col,1));
                        potential_synapses =  synapse.perms;
                        %ee = Ag(to_no,:)' == potential_synapses;
                        %aa = (potential_synapses == APstoredevo_prob(t_o,:));
                        active(:,c) = (potential_synapses(:,c) > connectedPerm);
                        %Overlap_c(:,c) = active(:,c).*c_mat;
                        %Overlap_c(:,c) = active(:,c).*c_mato;
                        Overlap_c(:,c) = c_mato.*active(:,c);
                        Overlap_cchek(:,c,t_o) = c_mato.*active(:,c); %Tensor
                        %Overlap_c(:,c) = Ag(to_no,:)' == potential_synapses;
                        Overlap_z = sum(Overlap_c);
                        %max_ov = roundn(max(Overlap_z),1);
                         
                    end

            neighbours_ov = Overlap_z;

            Overlap_z = Overlap_z.*max(boost);


        %% Inhibition:
            for c = 1:col_size
            [minLocalActivity_n,desired_localActivity_val] = kth_Score(neighbours_ov,...
            desired_localActivity);
                    if(Overlap_z(c) > stimulusThreshold && Overlap_z(c) > minLocalActivity_n)
                        %K_win_n = nonzeros(minLocalActivity_n);
                        activeColumns(c) = c;
                        activeColumns_z(c) = c;   

                    else

                        activeColumns(c) = 0;
                        activeColumns_z(c) = 0; 

                    end

            end

            activeColumns_n = nonzeros(activeColumns_z);
            i_sel = randperm(numel(activeColumns_n));


         %% Learning:

                col_win = Overlap_c(:,activeColumns_n);
                [ro_col_win,co_col_win] = size(col_win);

                col_win_perm = potential_synapses(:,activeColumns_n);

                active_col_win = active(:,activeColumns_n);

         %Update Permanences:
         for c = 1:co_col_win
           for s = 1:ro_col_win

            if(active_col_win(s,c) == 1) % --> Active

            synapse.perms(s,activeColumns_n(c)) = synapse.perms(s,activeColumns_n(c))+synPermActiveInc;
            synapse.perms(s,activeColumns_n(c)) = min(1.0,synapse.perms(s,activeColumns_n(c)));

            else

            synapse.perms(s,activeColumns_n(c)) = synapse.perms(s,activeColumns_n(c))-synPermInactiveDec;
            synapse.perms(s,activeColumns_n(c)) = max(1.0,synapse.perms(s,activeColumns_n(c)));

            end

            end

          end
                    activeColumns_list = activeColumns_z;
                    activeDutyCycle = updateActiveDutyCycle(activeColumns_list,it);
                    activeDutyCycleNeighbors = mean(activeDutyCycle);
                    boost = boostFunction(activeDutyCycle,activeDutyCycleNeighbors);
                    Overlap_z_n(it,:) = Overlap_z;
                    overlapDutyCycle= UpdateOverlapDutyCycle(Overlap_z_n, stimulusThreshold, it);
                    max_duty = maxActiveDutyCycle(activeDutyCycle);
                    minDutyCycle = 0.01*max_duty;

                    %logical-analysis of increasePermanences
                    overlapDutyCycle_bool = overlapDutyCycle < minDutyCycle;
                    increasePermanences = overlapDutyCycle_bool*0.1;
                    increasePermanences_n = increasePermanences == 0;
                    increasePermanences = increasePermanences + increasePermanences_n;

                    %Update Synapse.Permanences:
                    for cc = 1:ro_col
                       synapse.perms(cc,:) = synapse.perms(cc,:).*increasePermanences;
                    end

                 %Compute Receptive Field:   
                 inhibitionRadius = averageReceptiveFieldSize(synapse.perms,active);
                 inhibitionRadius = nan_filter(inhibitionRadius);
                 inhibitionRadius_mu = mean(inhibitionRadius);

                 desired_localActivity = roundn(inhibitionRadius_mu*col_size,0);

                 count = it;

        %% Extract all possible predictions(winner SDRs)
        Overlap_c_o = Overlap_c';
        %for ij = 1:numel(activeColumns_n)
            
        %Overlap_c_o_n(:,:,t_o) = Overlap_c_o; %Tensor
        Overlap_c_o_sp = Overlap_c_o(activeColumns_n,:);
    
        end
        %% Union Principle - First-Tier Prediction:
        [r_overlap, c_overlap] = size(Overlap_c_o_sp);
        if(t_o > 1 && t_o <= r_overlap)
            st1_union(t_o,:) =  or(Overlap_c_o_sp(t_o-1,:),...
                Overlap_c_o_sp(t_o,:));
        end
        sample_count = t_o   
        Agn_train(t_o,:) = Agn(t_o,:);   
    end
%% SPATIAL-POOLING ENDS HERE.
%%
%% FEEDFORWARD-TEMPORAL-CLASSIFICATION BEGINS HERE:    
        
    if(t_o > samples_n && t_o <=l_c_rangeo)
        %Overlap_c_o_sp_n(:,:,t_o) = Overlap_c_o_n(activeColumns_n,:,t_o);
        %koo = t_o;
        %Agn_test(t_o,:) = Agn(t_o,:);
        Agn_test= Agn((samples_n+1):l_c_rangeo,:);
        [rAgntest,cAgntest] = size(Agn_test);
        %[ro_sp,co_sp] = size(Overlap_c_o_sp);
        %Overlap_c_o_sp_nn = Overlap_c_o(activeColumns_n,:);
        [ro_sp,co_sp] = size(st1_union);
        
        for t_oo = 1:rAgntest
          
             %% Greedy Search:
             for i = 1: ro_sp

                 ee(i,1) = sum(Agn_test(t_oo,1:cAgntest-1)  == ...
                                st1_union(i,1:cAgntest-1));
                 if(ee(i,1) == max(ee))
                     sol_i(i,1) = i;
                 end

             end
             
             %% Classifier:
             sol_i = unique(nonzeros(sol_i));
             Overlap_c_o_sp_sol = st1_union(sol_i,:); 
             [ro_sp1,co_sp1] = size(Overlap_c_o_sp_sol);
             %ro_sp1;
              %% Union Principle - Second-Tier Prediction:
             for z = 2:ro_sp1
                  st_o_union(t_oo,:) = or(Overlap_c_o_sp_sol(z-1,:),...
                      Overlap_c_o_sp_sol(z,:));
             end
        end

    end
    
%% FEEDFORWARD-TEMPORAL-CLASSIFICATION ENDS HERE.
     
end
        

%% Compute Percentage Acurracy:
 class_acc = (sum(Agn_test(:,cAgntest)==st_o_union(:,cAgntest))...
     /(l_c_rangeo-samples_n))*100
 
 
%% Reference (BAMI book): 
%Hawkins, J., Ahmad, S., Purdy, S., & Lavin, A. (2016). 
%Biological and machine intelligence (Initial online release 0.4). Numenta Inc.
