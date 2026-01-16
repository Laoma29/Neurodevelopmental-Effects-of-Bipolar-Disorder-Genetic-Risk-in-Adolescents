%%双相情感障碍
clc;
clear  all;
close all;
% [data11,str11]=importdata('C:\Users\laoma\Desktop\合作课题\代码和数据\不对称性\DeviationZ_results\DeviationZ_results\allData\z-score_area_female.csv');
% [data22,str22]=importdata('C:\Users\laoma\Desktop\合作课题\代码和数据\不对称性\DeviationZ_results\DeviationZ_results\allData\z-score_area_male.csv');
[data_upps,str_upps]=importdata('C:\Users\laoma\Desktop\合作课题\代码和数据\不对称性\mh_y_pps.csv');

[data11,str11]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_thk_dsk.csv');
% [data22,str22]=importdata('C:\Users\laoma\Desktop\合作课题\代码和数据\不对称性\DeviationZ_results\DeviationZ_results\allData\z-score_area_male.csv');


textdata=data11.textdata(:,1);
cortex_thickness1=[];
score_upps1=[];
jk=1;
for i =1:length(data_upps.textdata)-1
text_str1=cell2mat(data_upps.textdata(i+1,1));
text_str2=cell2mat(data_upps.textdata(i+1,2));

for j=1:length(textdata)-1
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(text_str2,'baseline_year_1_arm_1')
Name_CT_Nor1(jk,:)=aa;
cortex_thickness1(jk,:)=data11.data(j,:);
jk=jk+1;

end
end
end

textdata=data11.textdata(:,1);
cortex_thickness2=[];
score_upps2=[];
jk=1;
for i =1:length(data_upps.textdata)-1
text_str1=cell2mat(data_upps.textdata(i+1,1));
text_str2=cell2mat(data_upps.textdata(i+1,2));
for j=1:length(textdata)-1
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(text_str2,'2_year_follow_up_y_arm_1')
Name_CT_Nor2(jk,:)=aa;
cortex_thickness2(jk,:)=data11.data(j,:);
jk=jk+1;

end
end
end



% [data33,str33]=importdata('ABCD_BP_PRSice_zscore.csv');
% 
% for i=1:length(data33.textdata)-1
% aaa=data33.textdata{1+i, 2};  
% PRS_zscore(i,:)=data33.data(i,2);
% PRS_raw(i,:)=data33.data(i,1);
% Name_RPS_group{i}=aaa(11:end);
% end
% 
% %%对齐数据
% kk=1;
% joint_var_matrix_final=[];
% score_upps_final=[];
% cortex_thickness_final=[];
% score_upps_z_final=[];
% PRS_zscore_final=[];
% PRS_raw_final=[];
% for i=1:length(Name_RPS_group)
% aa=cell2mat(Name_RPS_group(i));
% 
% for j=1:length(Name_CT_group)
% bb=Name_CT_group(j,:);
% 
% if strfind(aa,bb)
% Name_final(j,:)=bb;
% joint_var_matrix_final(:,:,kk)=joint_var_matrix(:,:,j);
% score_upps_final(kk,:)=score_upps(j,:);
% cortex_thickness_final(kk,:)=cortex_thickness(j,:);
% score_upps_z_final(kk,:)=score_upps_z(j,:);
% PRS_zscore_final(kk)=PRS_zscore(i,:);
% PRS_raw_final(kk)=PRS_raw(i,:);
% kk=kk+1;
% end
% end
% end





%%高低风险
load('Name_ID_Group.mat');
High_Risk_Name=Name_ID_Group.High_Risk_Name;
Low_Risk_Name=Name_ID_Group.Low_Risk_Name;
Risk_Name=[High_Risk_Name;Low_Risk_Name];

%%%高风险
textdata=data11.textdata(:,1);
textdata2=data11.textdata(:,2);
Name_CT_Nor_high_1=[];
cortex_thickness_high1=[];
score_upps1=[];
jk=1;
for i =1:length(High_Risk_Name)
text_str1=High_Risk_Name(i,:);
for j=1:length(textdata2)-1
    bb=cell2mat(textdata2(j+1,1));
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(bb,'baseline_year_1_arm_1')
Name_CT_Nor_high_1{jk}=aa;
cortex_thickness_high1(jk,:)=data11.data(j,:);
jk=jk+1;
end
end
end

Name_CT_Nor_high_2=[];
cortex_thickness_high2=[];
score_upps1=[];
jk=1;
for i =1:length(High_Risk_Name)
text_str1=High_Risk_Name(i,:);
for j=1:length(textdata2)-1
    bb=cell2mat(textdata2(j+1,1));
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(bb,'2_year_follow_up_y_arm_1')
Name_CT_Nor_high_2{jk}=aa;
cortex_thickness_high2(jk,:)=data11.data(j,:);
jk=jk+1;
end
end
end


%%%低风险
textdata=data11.textdata(:,1);
textdata2=data11.textdata(:,2);
Name_CT_Nor_low_1=[];
cortex_thickness_low1=[];
score_upps1=[];
jk=1;
for i =1:length(Low_Risk_Name)
text_str1=Low_Risk_Name(i,:);
for j=1:length(textdata2)-1
    bb=cell2mat(textdata2(j+1,1));
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(bb,'baseline_year_1_arm_1')
Name_CT_Nor_low_1{jk}=aa;
cortex_thickness_low1(jk,:)=data11.data(j,:);
jk=jk+1;
end
end
end


Name_CT_Nor_low_2=[];
cortex_thickness_low2=[];
score_upps1=[];
jk=1;
for i =1:length(Low_Risk_Name)
text_str1=Low_Risk_Name(i,:);
for j=1:length(textdata2)-1
    bb=cell2mat(textdata2(j+1,1));
    aa=cell2mat(textdata(j+1,1));
if strcmp(text_str1,aa) && strcmp(bb,'2_year_follow_up_y_arm_1')
Name_CT_Nor_low_2{jk}=aa;
cortex_thickness_low2(jk,:)=data11.data(j,:);
jk=jk+1;
end
end
end

%%%



%%%高风险
Name_CT_Nor_high1=[];
CT_all_high1=[];
CT_all_high2=[];
jk=1;
for i =1:length(Name_CT_Nor_high_1)
text_str1=cell2mat(Name_CT_Nor_high_1(i));
for j=1:length(Name_CT_Nor_high_2)
    aa=cell2mat(Name_CT_Nor_high_2(j));
if strcmp(text_str1,aa) 
Name_CT_Nor_high1(jk,:)=aa;
CT_all_high1(jk,:)=cortex_thickness_high1(i,:);
CT_all_high2(jk,:)=cortex_thickness_high2(j,:);
jk=jk+1;
end
end
end

%%%低风险
Name_CT_Nor_low1=[];
CT_all_low1=[];
CT_all_low2=[];
jk=1;
for i =1:length(Name_CT_Nor_low_1)
text_str1=cell2mat(Name_CT_Nor_low_1(i));
for j=1:length(Name_CT_Nor_low_2)
    aa=cell2mat(Name_CT_Nor_low_2(j));
if strcmp(text_str1,aa) 
Name_CT_Nor_low1(jk,:)=aa;
CT_all_low1(jk,:)=cortex_thickness_low1(i,:);
CT_all_low2(jk,:)=cortex_thickness_low2(j,:);
jk=jk+1;
end
end
end



CT_diff_group.Name_CT_Nor_low1=Name_CT_Nor_low1;
CT_diff_group.CT_all_low1=CT_all_low1;
CT_diff_group.CT_all_low2=CT_all_low2;
CT_diff_group.Name_CT_Nor_high1=Name_CT_Nor_high1;
CT_diff_group.CT_all_high1=CT_all_high1;
CT_diff_group.CT_all_high2=CT_all_high2;
%%
CT_all_high11=[];
CT_all_high22=[];
CT_all_low11=[];
CT_all_low22=[];


for k=1:length(CT_all_high1)
CT_all_high11(k,:)=zscore(CT_all_high1(k,1:68));
CT_all_high22(k,:)=zscore(CT_all_high2(k,1:68));
end

for k=1:length(CT_all_low1)
CT_all_low11(k,:)=zscore(CT_all_low1(k,1:68));
CT_all_low22(k,:)=zscore(CT_all_low2(k,1:68));
end

low_CT_diff=(CT_all_low22-CT_all_low11);%./CT_all_low11;%./cortex_thickness_low_final1;
high_CT_diff =(CT_all_high22-CT_all_high11);%./CT_all_high11;%./cortex_thickness_high_final1;

 [h,p,ci,stats] = ttest2(high_CT_diff,low_CT_diff);
% 
[inddd1,inddd2]=find(p<0.05);
stats.tstat(p<0.05)

b=low_CT_diff(:,inddd2(2))';%mean(Cn_Group(:,1:237),1)';
a=high_CT_diff(:,inddd2(2))';%mean(Cn_Group(:,238:end),1)';
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30];  
data2=[a b];
group_inx=[ones(length(a),1);zeros(length(b),1)];
h = daboxplot(data2,'groups',group_inx,...
    'colors',c,'whiskers',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'boxspacing',0.8); 




% Name_CT_Nor1;
% cortex_thickness1;
% 
% Name_CT_Nor2;
% cortex_thickness2;
% %%高风险
% %%第一次
% cortex_thickness_sub11=[];
% Name_high1=[];
% jk=1;
% for i =1:length(High_Risk_Name)
%     aa=High_Risk_Name(i,:);
%   for j=1:length(Name_CT_Nor1)  
%     bb=Name_CT_Nor1(j,:);
% if contains(aa,bb)
% Name_high1{jk}=aa;    
% cortex_thickness_sub11(jk,:)=cortex_thickness1(j,:);
% jk=jk+1;
% end
% end
% end
% 
% %%第二次
% Name_high_all=[];
% cortex_thickness_high11=[];
% cortex_thickness_high22=[];
% jk=1;
% for i =1:length(Name_high1)
%     aa=Name_high1{i};
%   for j=1:length(Name_CT_Nor2)  
%     bb=Name_CT_Nor2(j,:);
% if strcmp(aa,bb)
%  Name_high_all{jk} =bb; 
%  cortex_thickness_high11(jk,:)=cortex_thickness_sub11(i,:);
% cortex_thickness_high22(jk,:)=cortex_thickness2(j,:);
% jk=jk+1;
% end
% end
% end
% 
% %%低风险
% %%第一次
% cortex_thickness_sub11=[];
% Name_low1=[];
% jk=1;
% for i =1:length(Low_Risk_Name)
%     aa=Low_Risk_Name(i,:);
%   for j=1:length(Name_CT_Nor1)  
%     bb=Name_CT_Nor1(j,:);
% if contains(aa,bb)
% Name_low1{jk}=aa;    
% cortex_thickness_sub11(jk,:)=cortex_thickness1(j,:);
% jk=jk+1;
% end
% end
% end
% 
% %%第二次
% Name_low_all=[];
% cortex_thickness_low11=[];
% cortex_thickness_low22=[];
% jk=1;
% for i =1:length(Name_high1)
%     aa=Name_high1{i};
%   for j=1:length(Name_CT_Nor2)  
%     bb=Name_CT_Nor2(j,:);
% if strcmp(aa,bb)
%  Name_low_all{jk} =bb; 
%  cortex_thickness_low11(jk,:)=cortex_thickness_sub11(i,:);
% cortex_thickness_low22(jk,:)=cortex_thickness2(j,:);
% jk=jk+1;
% end
% end
% end
% 
% %%
% [~, idx_high] = unique(Name_high_all);  % 按行去重
% [~, idx_low] = unique(Name_low_all);  % 按行去重
% 
% %%
% cortex_thickness_low_final1=zscore(cortex_thickness_low11(idx_low,:));
% cortex_thickness_low_final2=zscore(cortex_thickness_low22(idx_low,:));
% %%
% cortex_thickness_high_final1=zscore(cortex_thickness_high11(idx_high,:));
% cortex_thickness_high_final2=zscore(cortex_thickness_high22(idx_high,:));
% 
% low_CT_diff=(cortex_thickness_low_final2-cortex_thickness_low_final1);%./cortex_thickness_low_final1;
% high_CT_diff =(cortex_thickness_high_final2-cortex_thickness_high_final1);%./cortex_thickness_high_final1;
% 
%  [h,p,ci,stats] = ttest2(high_CT_diff,low_CT_diff);
% 
% aaa=sum(score_upps_final(:,1:65),2);
% 
% bbb=mean(score_upps_z_final,2);
% 
% %%原始的结构数据
% % [data33,str33]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_thk_dsk.csv');
% % [data44,str44]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_vol_dsk.csv');
% % % [data44,str44]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_area_dsk.csv');
% % % [data55,str55]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_sulc_dsk.csv');
% % % [data66,str66]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_t1_gray_dsk.csv');
% % % [data77,str77]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\mri_y_smr_t1_white_dsk.csv');
% % % [data99,str99]=importdata('C:\Users\laoma\Desktop\合作课题\ABCD数据\smri\mri_y_smr_t1_gray_dsk.csv');
% % textdata=data33.textdata(:,1);
% % cortex_thickness11=[];
% % Cor_Vol11=[];
% % Name_CT_Group11=[];
% % jk=1;
% % for i =1:length(Name_final)
% % text_str1=cell2mat(data33.textdata(i+1,1));
% % text_str2=cell2mat(data33.textdata(i+1,2));
% % for j=1:length(textdata)-1
% %     aa=cell2mat(textdata(j+1,1));
% % if strcmp(text_str1,aa) && strcmp(text_str2,'baseline_year_1_arm_1')
% % Name_CT_Group11{jk}=aa;
% % cortex_thickness11(jk,:)=data33.data(j,:);
% % Cor_Vol11(jk,:)=data44.data(j,:);
% % 
% % jk=jk+1;
% % end
% % end
% % end
% % 
% % for l=1:size(joint_var_matrix,3)
% % cortex_thickness_new{l}=squeeze(joint_var_matrix(:,:,l)');
% % end
% % 
% % %%诊断
% % data_bip13 = readtable('C:\Users\laoma\Desktop\合作课题\代码和数据\不对称性\mh_y_ksads_bip.csv');
% % 
% % data_bip_cell13=table2cell(data_bip13);
% % data_bip_cell14=data_bip_cell13(:,[3:25 27:51 53 55 60:75 77 78 83 87:92]);
% % 
% % Data_bip_sub=zeros(size(data_bip_cell14));
% % for l=1:size(data_bip_cell14,1)
% %     for ll=1:size(data_bip_cell14,2)
% %         if isnan(data_bip_cell14{l,ll})
% %         Data_bip_sub(l,ll)=0;
% %         else
% %         Data_bip_sub(l,ll)=data_bip_cell14{l,ll};
% %         end
% %     end
% % end
% % 
% % %%对齐数据
% % kk=1;
% % joint_var_matrix_final13=[];
% % score_upps_final13=[];
% % cortex_thickness_final13=[];
% % score_upps_z_final13=[];
% % PRS_zscore_final13=[];
% % PRS_raw_final13=[];
% % 
% % for i=1:length(Name_final)
% % aa=(Name_final(i,:));
% % for j=1:length(Name_CT_group)
% % bb=Name_CT_group(j,:);
% % if strfind(aa,bb)
% % Name_final(j,:)=bb;
% % joint_var_matrix_final(:,:,kk)=joint_var_matrix(:,:,j);
% % score_upps_final(kk,:)=score_upps(j,:);
% % cortex_thickness_final(kk,:)=cortex_thickness(j,:);
% % score_upps_z_final(kk,:)=score_upps_z(j,:);
% % PRS_zscore_final(kk)=PRS_zscore(i,:);
% % PRS_raw_final(kk)=PRS_raw(i,:);
% % kk=kk+1;
% % end
% % end
% % end
% % 
% % Low_Risk_Name=Name_final(ind2(1:237),:);
% % High_Risk_Name=Name_final(ind2((end-237):end),:);
% % 
% % Name_ID_Group.Low_Risk_Name=Low_Risk_Name;
% % Name_ID_Group.High_Risk_Name=High_Risk_Name;
% % 
% % 
% % %%
% % score_upps13=[];
% % joint_var_matrix_final13=[];
% % score_upps_final13=[];
% % cortex_thickness_final13=[];
% % score_upps_z_final13=[];
% % PRS_zscore_final13=[];
% % PRS_raw_final13=[];
% % jkjk=1;
% % for j=1:length(Name_final)
% %     aa=(Name_final(j,:));
% % for i =1:length(data_bip_cell13)
% % text_str1=(data_bip_cell13{i,1});
% % text_str2=(data_bip_cell13{i,2});
% % data_bip_sub13=Data_bip_sub(i,:);
% % if strcmp(text_str1,aa) && strcmp(text_str2,'baseline_year_1_arm_1')
% % score_upps13(jkjk,:)=data_bip_sub13; 
% % joint_var_matrix_final13(:,:,jkjk)=joint_var_matrix_final(:,:,j);
% % score_upps_final13(jkjk,:)=score_upps_final(j,:);
% % cortex_thickness_final13(jkjk,:)=cortex_thickness_final(j,:);
% % score_upps_z_final13(jkjk,:)=score_upps_z_final(j,:);
% % PRS_zscore_final13(jkjk)=PRS_zscore_final(j);
% % PRS_raw_final13(jkjk)=PRS_raw_final(j);
% % jkjk=jkjk+1;
% % end
% % end
% % end
% % 
% % for i=1:size(joint_var_matrix_final13,3)
% %     W = threshold_proportional(joint_var_matrix_final13(:,:,i), 0.8);
% %     [Vulner_node] = gretna_vulnerability_weight(W);
% %     % [Vulner_node] = gretna_vulnerability(A);
% %     [averlocE, locEi] = gretna_node_local_efficiency_weight(W);
% %     [avercc, cci] = gretna_node_clustcoeff_weight(W,1);
% %     [averk, ki] = gretna_node_degree_weight(W);
% %     Vulner_node_group13(:,i)=Vulner_node;
% %     locEi_node_group13(:,i)=locEi';
% %     cci_node_group13(:,i)=cci';
% %     ki_node_group13(:,i)=ki';
% % end
% % 
% % %%
% % for jk=1:size(score_upps13,2)
% % score_upps13_zscore(:,jk)=zscore(score_upps13(:,jk));
% % end
% % 
% % [ind_Final1,ind_Final2]=sort(PRS_zscore_final13);
% % 
% % 
% % 
% % %% PRS_raw_final
% % [ind1,ind2]=sort(PRS_zscore_final);
% % number_sig=floor(0.05*length(PRS_zscore_final));
% % High_thre=ind1(end-number_sig);
% % Low_thre=ind1(number_sig);
% % % Mean_Risk=mean(PRS_raw_final);
% % % Sta_Risk=(std(PRS_raw_final));
% % Low_Risk_CT=cortex_thickness_final(1:237,:);
% % High_Risk_CT=cortex_thickness_final((end-237):end,:);
% % 
% % Low_Risk_joint_var_matrix=joint_var_matrix_final(:,:,1:237);
% % High_Risk_joint_var_matrix=joint_var_matrix_final(:,:,(end-237):end);
% % 
% % joint_var_matrix_1313=cat(3,Low_Risk_joint_var_matrix,High_Risk_joint_var_matrix);
% % 
% % for i=1:size(joint_var_matrix_1313,3)
% %     W = threshold_proportional(joint_var_matrix_1313(:,:,i), 0.8);
% %     [Vulner_node] = gretna_vulnerability_weight(W);
% %     Cn_Group(:,i)=clustering_coef_bu(W);
% % %     % [SWP,delta_C,delta_L] = small_world_propensity(W);
% % %     % W_Group(:,:,i)=W;
% % Vulner_node_Group(:,i)=Vulner_node;
% % [Len] = avg_path_matrix(W);
% % Len_Group(:,i)=Len;
% % end
% % 
% % [h,p11,ci,stats11] = ttest2(Len_Group(1:237),Len_Group(238:end));
% % 
% % [h,p111,ci,stats111] = ttest2(mean(Cn_Group(:,1:237),1)',mean(Cn_Group(:,238:end),1)');
% % 
% % b=Len_Group(1:237)';%mean(Cn_Group(:,1:237),1)';
% % a=Len_Group(238:end)';%mean(Cn_Group(:,238:end),1)';
% % c =  [0.45, 0.80, 0.69;...
% %       0.98, 0.40, 0.35;...
% %       0.55, 0.60, 0.79;...
% %       0.90, 0.70, 0.30];  
% % data2=[a;b];
% % group_inx=[ones(length(a),1);zeros(length(b),1)];
% % h = daboxplot(data2,'groups',group_inx,...
% %     'colors',c,'whiskers',0,...
% %     'scatter',1,'scattersize',15,'scatteralpha',0.5,...
% %     'boxspacing',0.8); 
% % 
% % 
% % 
% % 
% % for i=1:68
% %     [h,p1,ci,stats1] = ttest2(High_Risk_CT(:,i),Low_Risk_CT(:,i));
% % p_final(i,:)=p1;
% % Y_values(i)=stats1.tstat;
% % end
% %  [pID,pN] = FDR(p_final,0.05);
% % Y_values(p_final>pID)=0;
% % 
% % %%11 12
% % for i=1:68
% [ind1,ind2]=find(High_Risk_CT(:,i)>2);
% High_Posi_CT(i,:)=length(ind1)./number_sig;
% [ind1,ind2]=find(High_Risk_CT(:,i)<-2);
% High_Nega_CT(i,:)=length(ind1)./number_sig;
% 
% [ind1,ind2]=find(Low_Risk_CT(:,i)>2);
% Low_Posi_CT(i,:)=length(ind1)./number_sig;
% [ind1,ind2]=find(Low_Risk_CT(:,i)<-2);
% Low_Nega_CT(i,:)=length(ind1)./number_sig;
% end
% 
% 
% % Low_Risk_Neg_Num=sum(Low_Risk_CT<-2,1)./number_sig;
% % Low_Risk_Pos_Num=sum(Low_Risk_CT>2,1)./number_sig;
% % High_Risk_Neg_Num=sum(High_Risk_CT<-2,1)./number_sig;
% % High_Risk_Pos_Num=sum(High_Risk_CT>2,1)./number_sig;
% 
% %%
% 
% 
% % Low_Risk_Neg_CT=sum(Low_Risk_CT<-2,1)./number_sig;
% % Low_Risk_Pos_CT=sum(Low_Risk_CT>2,1)./number_sig;
% % High_Risk_Neg_CT=sum(High_Risk_CT<-2,1)./number_sig;
% % High_Risk_Pos_CT=sum(High_Risk_CT>2,1)./number_sig;
% 
% 
% p_final13=[];
% number_sig=910;
% YY_score_upps_z_final=mean(score_upps13_zscore,2);
% % [h,p22,ci,stats22] = ttest2(YY_score_upps_z_final(ind2(1:number_sig)),YY_score_upps_z_final(ind2((end-number_sig):end)));
% a=YY_score_upps_z_final(ind2(1:number_sig));
% b=YY_score_upps_z_final(ind2((end-number_sig):end));
% [h,p22,ci,stats22] = ttest2(a,b);
% 
% 
% c =  [0.45, 0.80, 0.69;...
%       0.98, 0.40, 0.35;...
%       0.55, 0.60, 0.79;...
%       0.90, 0.70, 0.30];  
% data2=[a;b];
% group_inx=[ones(length(a),1);zeros(length(b),1)];
% h = daboxplot(data2,'groups',group_inx,...
%     'colors',c,'whiskers',0,...
%     'scatter',1,'scattersize',15,'scatteralpha',0.5,...
%     'boxspacing',0.8); 
% 
% 
% for i=1:68
%     [h,p1,ci,stats1] = ttest2(cortex_thickness_final(ind2(1:number_sig),i),cortex_thickness_final(ind2((end-number_sig):end),i));
% p_final(i,:)=p1;
% Y_values(i)=stats1.tstat;
% end
%  [pID,pN] = FDR(p,0.05);
% Y_values(p>pID)=0;
% % [ind11,ind22]=sort(Y_values);
% % for jk=10:75
% % YY_score_upps_z_final=mean(score_upps13_zscore(:,1:jk),2);
% % [h,p22,ci,stats22] = ttest2(YY_score_upps_z_final(PRS_zscore_final13>=High_thre),YY_score_upps_z_final(PRS_zscore_final13<=Low_thre));
% % if p22<0.05
% % break;
% % end
% % end
% 
% 
% 
% shaffe_200=load('C:\Users\laoma\Desktop\Brainstorm_results\OMEGA_Dataset\SOurce_and_sinks\shaffe_200.mat');
% % % template=load('C:\Users\laoma\Desktop\Brainstorm_results\can_cam_dataset\FSAverage_2020\tess_cortex_pial_high.mat');
% [Newtemplate, I13, J] = MRI_tess_downsize_new( shaffe_200.template, 8001, 'reducepatch' );
% VertConn = tess_vertconn(Newtemplate.Vertices, Newtemplate.Faces);
% 
% [Vertices_sm, A] = tess_smooth(Newtemplate.Vertices, 0.15, 100, VertConn, 1, Newtemplate.Faces);
% Newtemplate.Vertices=Vertices_sm;
% [lH_model,rH_model]=split_cortical_model_new(Newtemplate);
% 
% 
% High_Posi_CT=sum(High_Risk_CT>2,1);
% High_Nega_CT=sum(High_Risk_CT<-2,1);
% Low_Posi_CT=sum(Low_Risk_CT>2,1);
% Low_Nega_CT=sum(Low_Risk_CT<-2,1);
% 
% 
% Data_final=zeros(8004,1);
% for i=1:68
% Data_final(Newtemplate.Atlas(3).Scouts(i).Vertices)=Low_Nega_CT(i);%stats1.tstat(i);
% end
% 
% h = view_surface_lxb('test',lH_model.Faces,lH_model.Vertices,Data_final(lH_model.lH),100);
% hold on;
% caxis([0 0.3]);
% % caxis([-3 3]);
% 
% hold on;
% color_number=100;
% CT=cbrewer('seq', 'BuPu', color_number);
% % CT=cbrewer('seq', 'YlOrBr', color_number);
% % CT=cbrewer('div', 'RdBu', color_number);
% CT(CT>1)=1;
% CT(CT<0)=0;
% colormap(CT);
% 
% h = view_surface_lxb('test',rH_model.Faces,rH_model.Vertices,Data_final(rH_model.rH),100);
% hold on;
% caxis([0 0.3]);
% % caxis([-3 3]);
% hold on;
% color_number=100;
% CT=cbrewer('div', 'RdBu', color_number);
% % CT=cbrewer('seq', 'BuPu', color_number);
% % CT=cbrewer('seq', 'YlOrBr', color_number);
% CT(CT>1)=1;
% CT(CT<0)=0;
% colormap(CT);
% % opts=struct('c',[],'maxiter',[],'PCAdim',[],'tol',1e-6,'epsilon',0.03);%% Change the epsilon value to .03- Suggested by Rajan for rsfMRI
% % [Ac1, Bc1] = cobe(FC_final,opts);
% 
% %%排序风险
% [ind1,ind2]=sort(PRS_zscore_final);
% number_par=floor(length(PRS_zscore_final)*0.05);
% case_pattern=cortex_thickness_final(1:number_par,:);
% disease_pattern=cortex_thickness_final((length(cortex_thickness_final)-number_par):end,:);
% 
% [h,p,ci,stats] = ttest2(disease_pattern,case_pattern);
% 
% [coef1, pval1] = corr(cortex_thickness_final, PRS_zscore_final');
% 
%     Vulner_node_group1=[];
%     locEi_node_group1=[];
%     cci_node_group1=[];
%     ki_node_group1=[];
% for i=1:size(joint_var_matrix_final13,3)
%     W = threshold_proportional(joint_var_matrix_final13(:,:,i), 0.8);
%     [Vulner_node] = gretna_vulnerability_weight(W);
%     % [Vulner_node] = gretna_vulnerability(A);
%     [averlocE, locEi] = gretna_node_local_efficiency_weight(W);
%     [avercc, cci] = gretna_node_clustcoeff_weight(W,1);
%     [averk, ki] = gretna_node_degree_weight(W);
%     Vulner_node_group1(:,i)=Vulner_node;
%     locEi_node_group1(:,i)=locEi';
%     cci_node_group1(:,i)=cci';
%     ki_node_group1(:,i)=ki';
% 
% 
% 
% %     W = threshold_proportional(abs(coef), 0.1);
% %     strength_Group(:,i)=strengths_und(W);
% %     degrees_Group(:,i)=degrees_und(W);
% %     Cn_Group(:,i)=clustering_coef_bu(W);
% %     % [SWP,delta_C,delta_L] = small_world_propensity(W);
% %     % W_Group(:,:,i)=W;
% %     %%
% % [Len] = avg_path_matrix(abs(coef));
% % 
% end
% 
% number_sig=237;
% [ind1,ind2]=sort(PRS_zscore_final13);
%  Vulner_node_group11=Vulner_node_group1';
% for i=1:68
% [h,p2,ci,stats2] = ttest2(Vulner_node_group11(ind2(1:number_sig),i),Vulner_node_group11(ind2(end-number_sig):end,i));
% p_final2(i,:)=p2;
% tstat_values(i,:)=stats2.tstat;
% end
% 
%  [pID2,pN] = FDR(p_final2,0.05);
% tstat_values(p_final2>pID2)=0;
% 
% Data_final=zeros(8004,1);
% for i=1:68
% Data_final(Newtemplate.Atlas(3).Scouts(i).Vertices)=tstat_values(i);%stats1.tstat(i);
% end
% h = view_surface_lxb('test',Newtemplate.Faces,Vertices_sm,Data_final,100);
% hold on;
% caxis([-3 3]);
% 
% % 
% % [h,p11,ci,stats11] = ttest2(disease_pattern,case_pattern);
% % 
% % [h,p22,ci,stats22] = ttest2(disease_pattern,case_pattern);
