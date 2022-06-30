clc
clear all
close all

%% Load sample data

load Dataset_sheet
load Dataset_signals

gestations_test = table2array(Dataset_sheet(:,1));
maternal_age = table2array(Dataset_sheet(:,2));
maternal_BMI = table2array(Dataset_sheet(:,4))./((table2array(Dataset_sheet(:,3))./100).^2);

%% Segmentation to one-minute segments
fs = 1000;
Coupling_labels = cell.empty;
Dataset_maternal_normal_perminute = cell.empty;
Dataset_fetal_normal_perminute = cell.empty;
Dataset_maternal_normal_beats_perminute = cell.empty;
Dataset_fetal_normal_beats_perminute = cell.empty;
for patient_id = 1:length(maternal_ECG)
    
    [patient_id length(maternal_ECG)]
        
selected_maternal_signal = maternal_ECG{patient_id,1};
selected_maternal_beats = maternal_beats{patient_id,1};

selected_fetal_signal = fetal_ECG{patient_id,1};
selected_fetal_beats = fetal_beats{patient_id,1};

%%%%%%%% to delete them in loop
maternal_beats_data = selected_maternal_beats;
fetal_beats_data = selected_fetal_beats;
time_length = length(selected_maternal_signal)/(1000*60); %minutes
one_minute = 60 * 1000;
full_minutes = 60 * time_length * 1000;
before = 1;
after = one_minute;
time = one_minute;
time_before = 0;
maternal_signal_oneminute = cell.empty;
fetal_signal_oneminute = cell.empty;
maternal_beats_oneminute = cell.empty;
fetal_beats_oneminute = cell.empty;
for j = 1:full_minutes/one_minute
        
    %%%%%%%%%%%%%%% Signals
    maternal_signal_oneminute{j,1} = selected_maternal_signal(before:after);
    fetal_signal_oneminute{j,1} = selected_fetal_signal(before:after);
    
    before = after + 1;
    after = after + one_minute;
    
    %%%%%%%%%%%%%%% Beats location
    location = find(maternal_beats_data <= time);
    location_last = location(end);
    location2 = find(fetal_beats_data <= time);
    location_last2 = location2(end);
    
    maternal_beats_oneminute{j,1} = maternal_beats_data(1:location_last) - time_before;
    fetal_beats_oneminute{j,1} = fetal_beats_data(1:location_last2) - time_before;
    
    if sum(isnan(maternal_beats_oneminute{j,1})) > 0
        [patient_id j]
    end
    if sum(isnan(fetal_beats_oneminute{j,1})) > 0
        [patient_id j]
    end
    
    maternal_beats_data(1:location_last) = [];
    fetal_beats_data(1:location_last2) = [];
    
    time_before = time;  
    time = time + one_minute;
    
end

overall_names = cell.empty;
for q = 1:time_length
    
selected_maternal_beats2 = maternal_beats_oneminute{q,1}';
selected_fetal_beats2 = fetal_beats_oneminute{q,1}';

frr = rand(400,1);
fs = 1000;

maternal_coupling_scenarios = 3;
names2 = ["",""];
for maternal_beats_number = 1:maternal_coupling_scenarios

[sF1,sM1,mphi,MtoF,fetal_beats_output,inst_fHR] = Phase_maternal_to_fetal(maternal_beats_number,selected_maternal_beats2',selected_fetal_beats2',frr');

fb_cnts=unique(fetal_beats_output);
focur = double.empty;
for j=1:length(fb_cnts)
    v1=fb_cnts(j);
    n=length(find(v1==fetal_beats_output));
    focur(j,:)=[v1 (n/length(fetal_beats_output))*100];
end

names = string;
for p = 1:size(focur,1)
names(p,1) = string([num2str(maternal_beats_number),'-',num2str(focur(p,1))]);
names(p,2) = string(num2str(focur(p,2)));
end
names2 = [names2;names];

end
names2(1,:) = [];

overall_names{q,1} = names2;

end

Coupling_labels{patient_id,1} = overall_names;

Dataset_maternal_normal_perminute{patient_id,1} = maternal_signal_oneminute;
Dataset_maternal_normal_beats_perminute{patient_id,1} = maternal_beats_oneminute;

Dataset_fetal_normal_perminute{patient_id,1} = fetal_signal_oneminute;
Dataset_fetal_normal_beats_perminute{patient_id,1} = fetal_beats_oneminute;

end

%% Prepare dataset with respect to coupling scenarios (selecting maximum)

final_maternal = cell.empty;
final_fetal = cell.empty;
final_maternal_beats = cell.empty;
final_fetal_beats = cell.empty;
final_labels = categorical;
overall_patients_count = double.empty;
overall_gestations = double.empty;
overall_maternal_age = double.empty;
overall_maternal_BMI = double.empty;
overall_health_fetal = double.empty;
overall_health_maternal = double.empty;
for pat_id = 1:length(Dataset_maternal_normal_perminute)
    
    [pat_id length(Dataset_maternal_normal_perminute)]

    selected_maternal = Dataset_maternal_normal_perminute{pat_id,1};
    selected_fetal = Dataset_fetal_normal_perminute{pat_id,1};
    selected_maternal_beats = Dataset_maternal_normal_beats_perminute{pat_id,1};
    selected_fetal_beats = Dataset_fetal_normal_beats_perminute{pat_id,1};
    
    selected_coupling = Coupling_labels{pat_id,1};
    
    overall_maternal = cell.empty;
    overall_fetal = cell.empty;
    overall_maternal_beats = cell.empty;
    overall_fetal_beats = cell.empty;
    overall_labels = categorical;
    for min_id = 1:length(selected_coupling)
        
        selected_minute_maternal = selected_maternal{min_id,1};
        selected_minute_fetal = selected_fetal{min_id,1};
    
        selected_minute_maternal_beats = selected_maternal_beats{min_id,1};
        selected_minute_fetal_beats = selected_fetal_beats{min_id,1};
        
        selected_minute_coupling = selected_coupling{min_id,1};
    
        coupling_percent = double(selected_minute_coupling(:,2));

        [value,id] = max(coupling_percent);
        doubled_maternal = selected_minute_maternal;
        doubled_fetal = selected_minute_fetal;
        doubled_maternal_beats = selected_minute_maternal_beats;
        doubled_fetal_beats = selected_minute_fetal_beats;
        
        doubled_label = categorical(selected_minute_coupling(id,1));
          
        overall_maternal = [overall_maternal;doubled_maternal];
        overall_fetal = [overall_fetal;doubled_fetal];
        overall_maternal_beats = [overall_maternal_beats;doubled_maternal_beats];
        overall_fetal_beats = [overall_fetal_beats;doubled_fetal_beats];
        
        overall_labels = [overall_labels;doubled_label];
        
    end
    
    final_maternal = [final_maternal;overall_maternal];
    final_fetal = [final_fetal;overall_fetal];    
    final_maternal_beats = [final_maternal_beats;overall_maternal_beats];
    final_fetal_beats = [final_fetal_beats;overall_fetal_beats];
    
    patients_count = repmat(pat_id,[length(overall_maternal_beats),1]);
    patient_gestation = repmat(gestations_test(pat_id),[length(overall_maternal_beats),1]);
    patient_maternal_age = repmat(maternal_age(pat_id),[length(overall_maternal_beats),1]);
    patient_maternal_BMI = repmat(maternal_BMI(pat_id),[length(overall_maternal_beats),1]);

    overall_patients_count = [overall_patients_count;patients_count];
    overall_gestations = [overall_gestations;patient_gestation];
    overall_maternal_age = [overall_maternal_age;patient_maternal_age];
    overall_maternal_BMI = [overall_maternal_BMI;patient_maternal_BMI];

    final_labels = [final_labels;overall_labels];
    
end
    
final_labels2 = reordercats(final_labels);

[final_labels3,order_id] = sort(final_labels2);

%% Preparation of arrays for deep learning
Final_dataset = cell.empty;
Final_dataset_beats = cell.empty;
for p = 1:length(final_maternal)
    
    [p length(final_maternal)]
    
    selected_maternal = final_maternal{p,1};
    selected_fetal = final_fetal{p,1};
    selected_maternal_beats = final_maternal_beats{p,1};
    selected_fetal_beats = final_fetal_beats{p,1};
    
    Final_dataset{p,1} = [selected_maternal;selected_fetal];
    Final_dataset_beats{p,1} = selected_maternal_beats;
    Final_dataset_beats{p,2} = selected_fetal_beats;

end

final_labels2_edit = double(final_labels2);
classes = unique(final_labels2_edit);
Final_dataset_edited = cell.empty;
Final_dataset_beats_edited = cell.empty;
Final_gestations = double.empty;
Final_maternal_age = double.empty;
Final_maternal_BMI = double.empty;
Final_patients_count = double.empty;
Final_labels = categorical;
for id = 1:length(classes)
    selected_class = classes(id,1);
    class_id = find(final_labels2_edit == selected_class);
    class_length = length(class_id);
    
    if class_length >= 8
       fixed_data = Final_dataset(class_id,1);
       fixed_beats = Final_dataset_beats(class_id,:);
       fixed_labels = final_labels2(class_id,1);
       fixed_patients_count = overall_patients_count(class_id,1);
       fixed_gestations_count = overall_gestations(class_id,1);
       fixed_maternal_age_count = overall_maternal_age(class_id,1);
       fixed_maternal_BMI_count = overall_maternal_BMI(class_id,1);

       Final_dataset_edited = [Final_dataset_edited;fixed_data];
       Final_dataset_beats_edited = [Final_dataset_beats_edited;fixed_beats];
       Final_patients_count = [Final_patients_count;fixed_patients_count];
       Final_labels = [Final_labels;fixed_labels];
       Final_gestations = [Final_gestations;fixed_gestations_count];
       Final_maternal_age = [Final_maternal_age;fixed_maternal_age_count];
       Final_maternal_BMI = [Final_maternal_BMI;fixed_maternal_BMI_count];
    end
    
end

Final_labels = removecats(Final_labels);

%% Deep learning prediction
X_test = Final_dataset_beats_edited;
X_test2 = single.empty;
for j = 1:size(X_test,1)
    selected_maternal_beats = X_test{j,1};
    vector = zeros(1,60000);
    vector(selected_maternal_beats) = 1;
    selected_fetal_beats = X_test{j,2};
    selected_fetal_beats(find(isnan(selected_fetal_beats) == 1)) = [];
    vector2 = zeros(1,60000);
    vector2(selected_fetal_beats) = 1;
    X_test2(:,:,:,j) = [single(vector);single(vector2)];
end
Y_test = Final_labels;

load trained_network

predictedLabels = classify(net,X_test2);
test_acc = (sum(categorical(predictedLabels) == Y_test)./numel(Y_test)).*100;

%% Geeting deep coherence and phase coherence
Final_patients_count_sorted = Final_patients_count;
Evaluation_set_sorted = X_test2;
Overall_nets_sorted2 = net;
Overall_predictions_sorted2 = predictedLabels;
Final_gestations_sorted = Final_gestations;

%%%% Deep coherence
patients = unique(Final_patients_count_sorted);
overall_RMSE = double.empty;
overall_heatmap_positive = zeros(1,5*60000);
sum_minutes_12 = zeros(5,60000);
sum_minutes_23 = zeros(5,60000);
sum_minutes_35 = zeros(5,60000);
counts_12 = zeros(5,1);
counts_23 = zeros(5,1);
counts_35 = zeros(5,1);
for pat_id = 1:length(patients)
    
    [pat_id length(patients)]
    selected_pat = patients(pat_id,1);
    selected_pat_location = find(Final_patients_count_sorted == selected_pat);
    selected_Evaluation_set_sorted = Evaluation_set_sorted(:,:,:,selected_pat_location);
    selected_Overall_nets_sorted2 = Overall_nets_sorted2;
    selected_Overall_predictions_sorted2 = Overall_predictions_sorted2(selected_pat_location,:);
    selected_true_label = Y_test(selected_pat_location,:);

    Final_gestations_sorted2 = Final_gestations_sorted(selected_pat_location,1);
    
    % Plot of Attention Vs Lambda
    len = length(selected_Overall_predictions_sorted2);
    overall_up = single.empty;
    sums = single.empty;

    for min_id = 1:len
        selected_net = selected_Overall_nets_sorted2;
        selected_beats = selected_Evaluation_set_sorted(:,:,:,min_id);
        selected_label = selected_Overall_predictions_sorted2(min_id,1);
        selected_true_label2 = selected_true_label(min_id,1);
        
        map = gradCAM(selected_net,selected_beats,selected_label,'FeatureLayer','relu_2');
        sums(min_id,1) = sum(sum(map));
        if sum(sum(map)) <= 100
        map = gradCAM(selected_net,selected_beats,selected_label,'FeatureLayer','relu_1');
        end
        
        env = single.empty;
        for o = 1:2
        [up,lo] = envelope(map(o,:),4000,'peak');
        env(o,:) = up;
        end
        env(:,1:700) = env(:,701:700+700);%map(:,701).*ones(2,700);
        env(:,end-700:end) = env(:,end-700-700:end-700);%map(:,end-701).*ones(2,701);

        map2 = sum(env)./2;   

        if (selected_label == '1-2')
           counts_12(min_id,1) = counts_12(min_id,1) + 1;
           sum_minutes_12(min_id,:) = sum_minutes_12(min_id,:) + map2;
        end
        
        if (selected_label == '2-3')
           counts_23(min_id,1) = counts_23(min_id,1) + 1;
           sum_minutes_23(min_id,:) = sum_minutes_23(min_id,:) + map2;
        end
        
        if (selected_label == '3-5')
           counts_35(min_id,1) = counts_35(min_id,1) + 1;
           sum_minutes_35(min_id,:) = sum_minutes_35(min_id,:) + map2;
        end

    end
    
end     
     
% Preparation and averaging
for o = 1:size(counts_35,1)
    number12 = counts_12(o,1);
    number23 = counts_23(o,1);
    number35 = counts_35(o,1);
    sum_minutes_12(o,:) = sum_minutes_12(o,:) ./ number12;
    sum_minutes_23(o,:) = sum_minutes_23(o,:) ./ number23;
    sum_minutes_35(o,:) = sum_minutes_35(o,:) ./ number35;
end

overall_sum_minutes = single.empty;
for g=1:5
    selected_sum = sum_minutes_12(g,:);
overall_sum_minutes = [overall_sum_minutes,selected_sum];
end
        more_id = find(overall_sum_minutes >= 1);
        for i = 1:length(more_id)
            selected_up_id = more_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.9 + (1-0.9).*rand(1,1);
        end
        less_id = find(overall_sum_minutes <= 0);
        for i = 1:length(less_id)
            selected_up_id = less_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.0 + (0.1-0.0).*rand(1,1);
        end
        up12 = movmean(overall_sum_minutes,30000,2);

overall_sum_minutes = single.empty;
for g=1:5
    selected_sum = sum_minutes_23(g,:);
overall_sum_minutes = [overall_sum_minutes,selected_sum];
end
        more_id = find(overall_sum_minutes >= 1);
        for i = 1:length(more_id)
            selected_up_id = more_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.9 + (1-0.9).*rand(1,1);
        end
        less_id = find(overall_sum_minutes <= 0);
        for i = 1:length(less_id)
            selected_up_id = less_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.0 + (0.1-0.0).*rand(1,1);
        end
        up23 = movmean(overall_sum_minutes,30000,2);

overall_sum_minutes = single.empty;
for g=1:5
    selected_sum = sum_minutes_35(g,:);
overall_sum_minutes = [overall_sum_minutes,selected_sum];
end
        more_id = find(overall_sum_minutes >= 1);
        for i = 1:length(more_id)
            selected_up_id = more_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.9 + (1-0.9).*rand(1,1);
        end
        less_id = find(overall_sum_minutes <= 0);
        for i = 1:length(less_id)
            selected_up_id = less_id(1,i);
            overall_sum_minutes(1,selected_up_id) = 0.0 + (0.1-0.0).*rand(1,1);
        end
        up35 = movmean(overall_sum_minutes,30000,2);    
        
%%%% Phase coherence 
sum_minutes_12 = zeros(5,60000);
sum_minutes_23 = zeros(5,60000);
sum_minutes_35 = zeros(5,60000);
counts_12 = zeros(5,1);
counts_23 = zeros(5,1);
counts_35 = zeros(5,1);      
for pat_id = 1:length(patients)
    
    [pat_id length(patients)]
    
    selected_pat = patients(pat_id,1);
    selected_pat_location = find(Final_patients_count_sorted == selected_pat);
    selected_Evaluation_set_sorted = Evaluation_set_sorted(:,:,:,selected_pat_location);
    selected_Overall_nets_sorted2 = Overall_nets_sorted2;
    selected_Overall_predictions_sorted2 = Overall_predictions_sorted2(selected_pat_location,:);
    selected_true_label = Y_test(selected_pat_location,:);

    Final_gestations_sorted2 = Final_gestations_sorted(selected_pat_location,1);
    
    scenarios = {'[1-1]','[1-2]','[1-3]','[1-4]',...
                 '[2-1]','[2-2]','[2-3]','[2-4]',...
                 '[3-1]','[3-2]','[3-3]','[3-4]','[3-5]'};
    selected_maternal_beats = Final_dataset_beats_edited(selected_pat_location,1);
    selected_fetal_beats = Final_dataset_beats_edited(selected_pat_location,2);
    len2 = length(selected_Overall_predictions_sorted2);
    overall_lam = single.empty;
    for min_id = 1:len2   
    maternal_beats = selected_maternal_beats{min_id,1};
    fetal_beats = selected_fetal_beats{min_id,1};
    selected_label = selected_Overall_predictions_sorted2(min_id,1);
    selected_true_label2 = selected_true_label(min_id,1);
    selected_maternal = maternal_beats./1000;
    selected_fetal = fetal_beats./1000;
    len2 = 60000;
    fs = 1000;
    t = (0:len2)/fs;
    N = 70/5; %Window for lambda calculation

    subi = 1;
    ratiotot = 5;
    clr = [238,34,104]./255;

all_phi = single.empty;
all_tim = single.empty;
all_lam = single.empty;
all_tlam = single.empty;   
count = 0;
for mb = 1:3 

    for i = 1:ratiotot 
        clear phiall phi tim lam tlam;
            
        fb = i+mb-1;

        % Relative phase (Psi) in the time window with respect to MECG
        [phiall,phi,tim] = phase_shift_A(selected_fetal,selected_maternal,fb,mb,t);
        
        % Phase coupling index (lambda):
        [lam] = lambda_A(phi,N); 
        lam_interp = resample(lam,60001,length(lam),0);
        lam_interp(end) = [];
        le = (length(tim)-length(lam)); 
        tlam = (tim(1:end-le));
        
        count = count + 1;
        all_lam(count,:) = lam_interp;
        
    end
              
end

if strcmp(string(selected_true_label2),["1-2"])
picked_lambda = all_lam(2,:);
counts_12(min_id,1) = counts_12(min_id,1) + 1;
sum_minutes_12(min_id,:) = sum_minutes_12(min_id,:) + picked_lambda;
end
if strcmp(string(selected_true_label2),["2-3"])
picked_lambda = all_lam(7,:);
counts_23(min_id,1) = counts_23(min_id,1) + 1;
sum_minutes_23(min_id,:) = sum_minutes_23(min_id,:) + picked_lambda;
end
if strcmp(string(selected_true_label2),["3-5"])
picked_lambda = all_lam(13,:);
counts_35(min_id,1) = counts_35(min_id,1) + 1;
sum_minutes_35(min_id,:) = sum_minutes_35(min_id,:) + picked_lambda;
end  
      
    end

end

% Preparation and averaging
for o = 1:size(counts_35,1)
    number12 = counts_12(o,1);
    number23 = counts_23(o,1);
    number35 = counts_35(o,1);
    sum_minutes_12(o,:) = sum_minutes_12(o,:) ./ number12;
    sum_minutes_23(o,:) = sum_minutes_23(o,:) ./ number23;
    sum_minutes_35(o,:) = sum_minutes_35(o,:) ./ number35;
end

overall_sum_minutes2 = single.empty;
for g=1:5
    selected_sum = sum_minutes_12(g,:);
overall_sum_minutes2 = [overall_sum_minutes2,selected_sum];
end
up122 = movmean(overall_sum_minutes2,30000,2);
        
overall_sum_minutes2 = single.empty;
for g=1:5
    selected_sum = sum_minutes_23(g,:);
overall_sum_minutes2 = [overall_sum_minutes2,selected_sum];
end
up232 = movmean(overall_sum_minutes2,30000,2);
        
overall_sum_minutes2 = single.empty;
for g=1:5
    selected_sum = sum_minutes_35(g,:);
overall_sum_minutes2 = [overall_sum_minutes2,selected_sum];
end
up352 = movmean(overall_sum_minutes2,30000,2);
     
        len = 5;
        f = figure;%('Position',[443,42,866,954]);
        hold on;
        subplot(3,3,1),plot(up12,'color',[0.64,0.08,0.18],'LineWidth',2)
        hold on
        plot(up122,'--','color','r','LineWidth',2)
        rmse = sqrt(nanmean((up12 - up122).^2));
        scatter(-1000000,rmse,'filled','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[0 0 0],'LineWidth',2);
        xlim([0 len*60000]);
        xticks([0:floor(len*60000/len):len*60000]);
        tick_names = cell.empty;
        tick_names{1,1} = string(0);
        for tick_id = 1:len
            tick_names{1,tick_id+1} = string(tick_id);
        end
        xticklabels(tick_names);
        xlabel('Time (min)','FontWeight','bold');
        ylim([0 1.01]);
        yticks([0:0.2:1]);
        ylabel('Coupling strength','FontWeight','bold')
        grid on
        box on
        ax = gca;
        ax.FontSize = 10;
        xtickangle(0)
        lgd = legend('DL attention','Coupling Lambda',['RMSE: ',num2str(round(rmse,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',9,...
        'FontWeight','bold');
   
        subplot(3,3,4),plot(up23,'color',[0.00,0.42,0.00],'LineWidth',2)
        hold on
        plot(up232,'--','color','g','LineWidth',2)
        rmse = sqrt(nanmean((up23 - up232).^2));
        scatter(-1000000,rmse,'filled','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[0 0 0],'LineWidth',2);
        xlim([0 len*60000]);
        xticks([0:floor(len*60000/len):len*60000]);
        tick_names = cell.empty;
        tick_names{1,1} = string(0);
        for tick_id = 1:len
            tick_names{1,tick_id+1} = string(tick_id);
        end
        xticklabels(tick_names);
        xlabel('Time (min)','FontWeight','bold');
        ylim([0 1.01]);
        yticks([0:0.2:1]);
        ylabel('Coupling strength','FontWeight','bold')
        grid on
        box on
        ax = gca;
        ax.FontSize = 10;
        xtickangle(0)
        lgd = legend('DL attention','Coupling Lambda',['RMSE: ',num2str(round(rmse,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',9,...
        'FontWeight','bold');
    
        subplot(3,3,7),plot(up35,'color','b','LineWidth',2)
        hold on;
        plot(up352,'--','color',[0.07,0.62,1.00],'LineWidth',2) 
        rmse = sqrt(nanmean((up35 - up352).^2));
        scatter(-1000000,rmse,'filled','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[0 0 0],'LineWidth',2);
        xlim([0 len*60000]);
        xticks([0:floor(len*60000/len):len*60000]);
        tick_names = cell.empty;
        tick_names{1,1} = string(0);
        for tick_id = 1:len
            tick_names{1,tick_id+1} = string(tick_id);
        end
        xticklabels(tick_names);
        xlabel('Time (min)','FontWeight','bold');
        ylim([0 1.01]);
        yticks([0:0.2:1]);
        ylabel('Coupling strength','FontWeight','bold')
        grid on
        box on
        ax = gca;
        ax.FontSize = 10;
        xtickangle(0)
        lgd = legend('DL attention','Coupling Lambda',['RMSE: ',num2str(round(rmse,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',9,...
        'FontWeight','bold');
        
subplot(3,3,2),image([1 len*60000],[-5 7],up12,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)

subplot(3,3,5),image([1 len*60000],[-5 7],up23,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)

subplot(3,3,8),image([1 len*60000],[-5 7],up35,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)


subplot(3,3,3),image([1 len*60000],[-5 7],up122,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)

subplot(3,3,6),image([1 len*60000],[-5 7],up232,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)

subplot(3,3,9),image([1 len*60000],[-5 7],up352,'CDataMapping','scaled',...
                     'AlphaData',1);
xlim([0 len*60000]);
xticks([0:floor(len*60000/len):len*60000]);
tick_names = cell.empty;
tick_names{1,tick_id} = string(0);
for tick_id = 1:len
    tick_names{1,tick_id+1} = string(tick_id);
end
xticklabels(tick_names);
xlabel('Time (min)','FontWeight','bold');
ylim([0 1.01]);
box on
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 10;
caxis([0 1]);
colormap jet
c = colorbar('FontSize',10,'FontWeight','bold','Location','northoutside',...
    'Ticks',[0:0.2:1]);
c.Label.String = 'Coupling strength';
c.Label.FontSize = 10;
xtickangle(0)

% exportgraphics(f,'overall_heatmap_test.png');

%% Bland-Altman analysis
up12_2 = up12;
up122_2 = up122;
up12_edit = resample(double(up12_2),250,5*60000,0);
up122_edit = resample(double(up122_2),250,5*60000,0);

up23_2 = up23;
up232_2 = up232;
up23_edit = resample(double(up23_2),250,5*60000,0);
up232_edit = resample(double(up232_2),250,5*60000,0);

up35_edit = resample(double(up35),250,5*60000,0);
up352_edit = resample(double(up352),250,5*60000,0);
nan_index = find(isnan(up352_edit) == 1);
nan_index = nan_index(end);
up35_edit = up35_edit(nan_index+1:end);
up352_edit = up352_edit(nan_index+1:end);


fig = figure('Position',[463,323,887,493]);
ax = gca;
ax.FontSize = 20;
hold;
xtickangle(0);
ytickangle(0);
xlabel('Average','FontWeight','bold')
ylabel('Difference','FontWeight','bold')
grid on
box on
[means,diffs,meanDiff,CR,linFit] = BlandAltman(up12_edit, up122_edit, 3);
hold on
plot(means,diffs,'.k','MarkerSize',7)
plot(means, ones(1,length(means)).*CR(1),'-k','LineWidth',2); %%%plot the upper CR
plot(means, ones(1,length(means)).*CR(2),'-k','LineWidth',2); %%%plot the lower CR
plot(means,meanDiff.*ones(1,length(means)),'-k','LineWidth',2); %%%plot zero
minimum = min(means);
maximum = max(means);
plot([minimum,maximum], [minimum,maximum].*linFit(1)+linFit(2),'-r','LineWidth',2);
y1 = ones(1,length(means)).*CR(1);
text(nanmin(means),... 
     y1(1)+0.008,...
     'Mean+2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
y2 = ones(1,length(means)).*CR(2);
text(nanmin(means),... 
     y2(1)-0.05,...
     'Mean-2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
text(nanmin(means),... 
     meanDiff+0.008,...
     {['Mean= ',num2str(round(meanDiff,2))],['2Std= \pm',num2str(round(CR(1),2))]},...
     'HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');

fig = figure('Position',[463,323,887,493]);
ax = gca;
ax.FontSize = 20;
hold;
xtickangle(0);
ytickangle(0);
xlabel('Average','FontWeight','bold')
ylabel('Difference','FontWeight','bold')
grid on
box on
[means,diffs,meanDiff,CR,linFit] = BlandAltman(up23_edit, up232_edit, 3);
hold on
plot(means,diffs,'.k','MarkerSize',7)
plot(means, ones(1,length(means)).*CR(1),'-k','LineWidth',2); %%%plot the upper CR
plot(means, ones(1,length(means)).*CR(2),'-k','LineWidth',2); %%%plot the lower CR
plot(means,meanDiff.*ones(1,length(means)),'-k','LineWidth',2); %%%plot zero
minimum = min(means);
maximum = max(means);
plot([minimum,maximum], [minimum,maximum].*linFit(1)+linFit(2),'-r','LineWidth',2);
y1 = ones(1,length(means)).*CR(1);
text(nanmin(means),... 
     y1(1)+0.008,...
     'Mean+2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
y2 = ones(1,length(means)).*CR(2);
text(nanmin(means),... 
     y2(1)-0.05,...
     'Mean-2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
text(nanmin(means),... 
     meanDiff+0.008,...
     {['Mean= ',num2str(round(meanDiff,2))],['2Std= \pm',num2str(round(CR(1),2))]},...
     'HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');

fig = figure('Position',[463,323,887,493]);
ax = gca;
ax.FontSize = 20;
hold;
xtickangle(0);
ytickangle(0);
xlabel('Average','FontWeight','bold')
ylabel('Difference','FontWeight','bold')
grid on
box on
[means,diffs,meanDiff,CR,linFit] = BlandAltman(up35_edit, up352_edit, 3);
hold on
plot(means,diffs,'.k','MarkerSize',7)
plot(means, ones(1,length(means)).*CR(1),'-k','LineWidth',2); %%%plot the upper CR
plot(means, ones(1,length(means)).*CR(2),'-k','LineWidth',2); %%%plot the lower CR
plot(means,meanDiff.*ones(1,length(means)),'-k','LineWidth',2); %%%plot zero
minimum = min(means);
maximum = max(means);
plot([minimum,maximum], [minimum,maximum].*linFit(1)+linFit(2),'-r','LineWidth',2);
y1 = ones(1,length(means)).*CR(1);
text(nanmin(means),... 
     y1(1)+0.008,...
     'Mean+2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
y2 = ones(1,length(means)).*CR(2);
text(nanmin(means),... 
     y2(1)-0.05,...
     'Mean-2Std','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');
text(nanmin(means),... 
     meanDiff+0.008,...
     {['Mean= ',num2str(round(meanDiff,2))],['2Std= \pm',num2str(round(CR(1),2))]},...
     'HorizontalAlignment','left',...
     'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold',...
     'BackgroundColor','w');

%% With gestation

%%%%%%% Deep Coherence
patients = unique(Final_patients_count_sorted);
overall_RMSE = double.empty;
overall_heatmap_positive = zeros(1,600000);
sum_minutes_12 = nan(5,60000);
sum_minutes_23 = nan(5,60000);
sum_minutes_35 = nan(5,60000);
sum_gestations_12 = nan(5,60000);
sum_gestations_23 = nan(5,60000);
sum_gestations_35 = nan(5,60000);
counts_12 = zeros(5,1);
counts_23 = zeros(5,1);
counts_35 = zeros(5,1);
probabilitites_patients = double.empty;
for pat_id = 1:length(patients)
    
    [pat_id length(patients)]
    selected_pat = patients(pat_id,1);
    selected_pat_location = find(Final_patients_count_sorted == selected_pat);
    selected_Evaluation_set_sorted = Evaluation_set_sorted(:,:,:,selected_pat_location);
    selected_Overall_nets_sorted2 = Overall_nets_sorted2;
    selected_Overall_predictions_sorted2 = Overall_predictions_sorted2(selected_pat_location,:);
    selected_true_label = Y_test(selected_pat_location,:);

    Final_gestations_sorted2 = Final_gestations_sorted(selected_pat_location,1);

    prob_12 = length(find(selected_Overall_predictions_sorted2 == '1-2'))/length(selected_Overall_predictions_sorted2).*100;
    prob_23 = length(find(selected_Overall_predictions_sorted2 == '2-3'))/length(selected_Overall_predictions_sorted2).*100;
    prob_35 = length(find(selected_Overall_predictions_sorted2 == '3-5'))/length(selected_Overall_predictions_sorted2).*100;

    % Plot of Attention Vs Lambda
    len = length(selected_Overall_predictions_sorted2);
    overall_up = single.empty;
    sums = single.empty;

    for min_id = 1:len

        selected_net = selected_Overall_nets_sorted2;
        selected_beats = selected_Evaluation_set_sorted(:,:,:,min_id);
        selected_label = selected_Overall_predictions_sorted2(min_id,1);
        selected_true_label2 = selected_true_label(min_id,1);
        gestation_minutes = Final_gestations_sorted2(min_id,1);
        
        map = gradCAM(selected_net,selected_beats,selected_label,'FeatureLayer','relu_2');
        sums(min_id,1) = sum(sum(map));
        if sum(sum(map)) <= 100
        map = gradCAM(selected_net,selected_beats,selected_label,'FeatureLayer','relu_1');
        end
        
        env = single.empty;
        for o = 1:2
        [up,lo] = envelope(map(o,:),4000,'peak');
        env(o,:) = up;
        end
        env(:,1:700) = env(:,701:700+700);%map(:,701).*ones(2,700);
        env(:,end-700:end) = env(:,end-700-700:end-700);%map(:,end-701).*ones(2,701);

        map2 = sum(env)./2;   
        
        map2_average = mean(map2);
        
        if (selected_label == '1-2')
           counts_12(min_id,1) = counts_12(min_id,1) + 1;
           sum_minutes_12(min_id,counts_12(min_id,1)) = map2_average;
           sum_gestations_12(min_id,counts_12(min_id,1)) = gestation_minutes;
        end
        
        if (selected_label == '2-3')
           counts_23(min_id,1) = counts_23(min_id,1) + 1;
           sum_minutes_23(min_id,counts_23(min_id,1)) = map2_average;
           sum_gestations_23(min_id,counts_23(min_id,1)) = gestation_minutes;
        end
        
        if (selected_label == '3-5')
           counts_35(min_id,1) = counts_35(min_id,1) + 1;
           sum_minutes_35(min_id,counts_35(min_id,1)) = map2_average;
           sum_gestations_35(min_id,counts_35(min_id,1)) = gestation_minutes;
        end

    end
    
    probabilitites_patients(pat_id,:) = [prob_12;prob_23;prob_35];
    
end     

figure('Position',[463,323,887,493]);
hold on
gest_scenarios_12 = unique(sum_gestations_12);
gest_scenarios_12(isnan(gest_scenarios_12) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_12)
    selected_gest_scenario = gest_scenarios_12(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_12 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_12(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',clr,'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Deep learning','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);

figure('Position',[463,323,887,493]);
hold on
gest_scenarios_23 = unique(sum_gestations_23);
gest_scenarios_23(isnan(gest_scenarios_23) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_23)
    selected_gest_scenario = gest_scenarios_23(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_23 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_23(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',clr,'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Deep learning','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);
       
figure('Position',[463,323,887,493]);
hold on
gest_scenarios_35 = unique(sum_gestations_35);
gest_scenarios_35(isnan(gest_scenarios_35) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_35)
    selected_gest_scenario = gest_scenarios_35(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_35 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_35(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',clr,'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Deep learning','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);


%%%%%%% Phase Coherence
patients = unique(Final_patients_count_sorted);
overall_RMSE = double.empty;
overall_heatmap_positive = zeros(1,600000);
sum_minutes_12 = nan(5,60000);
sum_minutes_23 = nan(5,60000);
sum_minutes_35 = nan(5,60000);
sum_gestations_12 = nan(5,60000);
sum_gestations_23 = nan(5,60000);
sum_gestations_35 = nan(5,60000);
counts_12 = zeros(5,1);
counts_23 = zeros(5,1);
counts_35 = zeros(5,1);
probabilitites_patients = double.empty;
for pat_id = 1:length(patients)
    
    [pat_id length(patients)]
    selected_pat = patients(pat_id,1);
    selected_pat_location = find(Final_patients_count_sorted == selected_pat);
    selected_Evaluation_set_sorted = Evaluation_set_sorted(:,:,:,selected_pat_location);
    selected_Overall_nets_sorted2 = Overall_nets_sorted2;
    selected_Overall_predictions_sorted2 = Overall_predictions_sorted2(selected_pat_location,:);
    selected_true_label = Y_test(selected_pat_location,:);

    Final_gestations_sorted2 = Final_gestations_sorted(selected_pat_location,1);

    prob_12 = length(find(selected_Overall_predictions_sorted2 == '1-2'))/length(selected_Overall_predictions_sorted2).*100;
    prob_23 = length(find(selected_Overall_predictions_sorted2 == '2-3'))/length(selected_Overall_predictions_sorted2).*100;
    prob_35 = length(find(selected_Overall_predictions_sorted2 == '3-5'))/length(selected_Overall_predictions_sorted2).*100;

    scenarios = {'[1-1]','[1-2]','[1-3]','[1-4]',...
                 '[2-1]','[2-2]','[2-3]','[2-4]',...
                 '[3-1]','[3-2]','[3-3]','[3-4]','[3-5]'};
    selected_maternal_beats = Final_dataset_beats_edited(selected_pat_location,1);
    selected_fetal_beats = Final_dataset_beats_edited(selected_pat_location,2);
    len2 = length(selected_Overall_predictions_sorted2);
    overall_lam = single.empty;
    for min_id = 1:len2   
    maternal_beats = selected_maternal_beats{min_id,1};
    fetal_beats = selected_fetal_beats{min_id,1};
    selected_label = selected_Overall_predictions_sorted2(min_id,1);
    selected_true_label2 = selected_true_label(min_id,1);
    selected_maternal = maternal_beats./1000;
    selected_fetal = fetal_beats./1000;
    gestation_minutes = Final_gestations_sorted2(min_id,1);
    
    len2 = 60000;
    fs = 1000;
    t = (0:len2)/fs;
    N = 70/5;

    subi = 1;
    ratiotot = 5;
    clr = [238,34,104]./255;

all_phi = single.empty;
all_tim = single.empty;
all_lam = single.empty;
all_tlam = single.empty;   
count = 0;
for mb = 1:3 

    for i = 1:ratiotot 
        clear phiall phi tim lam tlam;
            
        fb = i+mb-1;

        % Relative phase (Psi) in the time window with respect to MECG
        [phiall,phi,tim] = phase_shift_A(selected_fetal,selected_maternal,fb,mb,t);
        
        % Phase coupling index (lambda):
        [lam] = lambda_A(phi,N); 
        lam_interp = resample(lam,91,length(lam),0);
        lam_interp(end) = [];
        le = (length(tim)-length(lam)); 
        tlam = (tim(1:end-le));
        
        count = count + 1;
        all_lam(count,:) = lam_interp;
        
    end
              
end

all_lam(all_lam >= 0.95) = 0.9;
all_lam(all_lam <= 0.2) = 0.2;
if strcmp(string(selected_true_label2),["1-2"])
picked_lambda = all_lam(2,:);
picked_lambda = movmean(picked_lambda,91);
picked_lambda_value = mean(picked_lambda);
% picked_lambda_value = quantile(all_lam(2,:),0.75);
counts_12(min_id,1) = counts_12(min_id,1) + 1;
sum_minutes_12(min_id,counts_12(min_id,1)) = picked_lambda_value;
sum_gestations_12(min_id,counts_12(min_id,1)) = gestation_minutes;
end
if strcmp(string(selected_true_label2),["2-3"])
picked_lambda = all_lam(7,:);
picked_lambda = movmean(picked_lambda,91);
picked_lambda_value = mean(picked_lambda);
% picked_lambda_value = quantile(all_lam(7,:),0.75);
counts_23(min_id,1) = counts_23(min_id,1) + 1;
sum_minutes_23(min_id,counts_23(min_id,1)) = picked_lambda_value;
sum_gestations_23(min_id,counts_23(min_id,1)) = gestation_minutes;
end
if strcmp(string(selected_true_label2),["3-5"])
picked_lambda = all_lam(13,:);
picked_lambda = movmean(picked_lambda,91);
picked_lambda_value = mean(picked_lambda);
% picked_lambda_value = quantile(all_lam(13,:),0.75);
counts_35(min_id,1) = counts_35(min_id,1) + 1;
sum_minutes_35(min_id,counts_35(min_id,1)) = picked_lambda_value;
sum_gestations_35(min_id,counts_35(min_id,1)) = gestation_minutes;
end  
      
    end
    
    probabilitites_patients(pat_id,:) = [prob_12;prob_23;prob_35];

end     

figure('Position',[463,323,887,493]);
hold on
gest_scenarios_12 = unique(sum_gestations_12);
gest_scenarios_12(isnan(gest_scenarios_12) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_12)
    selected_gest_scenario = gest_scenarios_12(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_12 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_12(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',[0.07,0.62,1.00],'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Lambda','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);

figure('Position',[463,323,887,493]);
hold on
gest_scenarios_23 = unique(sum_gestations_23);
gest_scenarios_23(isnan(gest_scenarios_23) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_23)
    selected_gest_scenario = gest_scenarios_23(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_23 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_23(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',[0.07,0.62,1.00],'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Lambda','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);

figure('Position',[463,323,887,493]);
hold on
gest_scenarios_35 = unique(sum_gestations_35);
gest_scenarios_35(isnan(gest_scenarios_35) == 1) = [];
vector_gestations = single.empty;
vector_attentions = single.empty;
for gest_id = 1:length(gest_scenarios_35)
    selected_gest_scenario = gest_scenarios_35(gest_id,1);
    [gest_row,gest_col] = find(sum_gestations_35 == selected_gest_scenario);
    selected_attentions = single.empty;
    for i = 1:length(gest_row)
        selected_attentions(i,1) = sum_minutes_35(gest_row(i,1),gest_col(i,1));
    end
    vector_gestations = [vector_gestations,repmat(selected_gest_scenario,1,length(selected_attentions))];
    vector_attentions = [vector_attentions,selected_attentions'];
end
vector_attentions(vector_attentions > 1) = 0.9;
vector_attentions(vector_attentions <= 0) = mean(vector_attentions);
sdDiff = nanstd(vector_attentions);
plot(vector_gestations,vector_attentions,'LineStyle','none','MarkerSize',8,'Marker','square',...
     'MarkerFaceColor',[0.07,0.62,1.00],'MarkerEdgeColor','k');
[linFit,S] = polyfit(vector_gestations,vector_attentions,1);
[y_fit,delta] = polyval(linFit,vector_gestations,S);
plot(vector_gestations, y_fit,'-k','LineWidth',2);
plot(vector_gestations, y_fit+2*delta,'--k','LineWidth',2); %%%plot the upper CR
plot(0,10,'.k','MarkerSize',20);
[rmcc,pmcc]=corrcoef(vector_gestations,vector_attentions); %m for matrix
R = rmcc(1,2);
p = pmcc(1,2);
plot(0,10,'.k','MarkerSize',20);
plot(0,10,'.k','MarkerSize',20);
plot(vector_gestations, y_fit-2*delta,'--k','LineWidth',2); %%%plot the lower CR
xlim([19 41]);
xticks([20:2:40]);
xlabel('Gestational age (weeks)','FontWeight','bold');
ylim([0 1.01]);
yticks([0:0.2:1]);
ylabel('Coupling','FontWeight','bold')
box on
grid on
ax = gca;
ax.FontSize = 20;
xtickangle(0)
lgd = legend('Lambda','Linear fitting','95% CI',...
             ['y = ',num2str(round(linFit(1,1),3)),'x + ',num2str(round(linFit(1,2),3))],...
             ['R = ',num2str(round(R,3))],...
             ['{\it p}-value = ',num2str(round(p,3))],...
        'Location','northoutside','Orientation','horizontal','FontSize',18,...
        'FontWeight','bold','NumColumns',3);



