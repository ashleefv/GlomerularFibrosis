
function gluc = GlucoseInputFunc(t,GlucoseCtrlOptn,TreatmentTIme)

    mol_2_gram = 180/1000000;                                    % Converting glucose units from mmol/l to g/ml
    weeks_2_days = 7;                                            % Converting the week inputs into days 
    y1 = 5.8;
    y2 = 28.3;
    t1 = 6*weeks_2_days;
    t2 = 16*weeks_2_days;
    TreatmentTimeDays = TreatmentTIme*weeks_2_days;

    if t <= t1
        gluc = y1*mol_2_gram;
    elseif t > t1 && t <= t2 
        gluc = (((y2-y1)/(t2-t1))*t + y1 - ((y2-y1)/(t2-t1))*t1)*mol_2_gram;
    elseif t > t2 && t <= TreatmentTimeDays
        gluc = y2*mol_2_gram;
    elseif t > TreatmentTimeDays
        if sum(strcmp(GlucoseCtrlOptn,'NoGlucoseCtrl'))>0
            gluc = y2*mol_2_gram;
        elseif sum(strcmp(GlucoseCtrlOptn,'YesGlucoseCtrl'))>0
            gluc = y1*mol_2_gram;
        end
    end

end