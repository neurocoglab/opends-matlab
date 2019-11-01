idx_msg_condit = false(length(Time),1);

TF = TimeMsg(idx_condit);

j = 1;
for i = 1 : length(TF)
    
    tm = TF(i);
    
    while j <= length(Time) && tm > Time(j)
       j = j + 1;
    end
    
    if j <= length(Time)
       idx_msg_condit(j)=1; 
    end
    
end