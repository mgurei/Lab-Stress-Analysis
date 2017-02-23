%------------------QRS Detection Function---------------------------------%
%- This is a starting point including the filtering stages for a QRS
%algorithm. The b and a parameters are meant for 4 filtering stages, and
%can be found using the fdatool in Matlab. 

% INPUT
% --- ecg - ECG signal for QRS detection
% --- b_low - FIR low pass b filtering coefficients
% --- b_high - FIR high pass b filtering coefficients
% --- b_avg - averageing filter coefficients
% --- delay - delay caused by filters (check with grpdelay function)
function [qrs]=s_detect_qrs(ecg,b_low,b_high,b_avg,delay)
    
    % Important Values
    % Try to make improvals based on these variables or improve the
    % detection process itself
    fs=256;
    T_ref=0.25*fs; % 0.25 second refractory period
    window=3*fs; % 3 second window
    h_thresh=0; % initial value of h_thresh
    h_thresh_correct=0.4; % correction value for h_thresh
    
    % Detecting candidate 
    candidate_detected=0;
    candidate_pos=0;
    candidate=0;
    
    % Keeping track of maximum values in last 5 windows
    window_max_buff=zeros(1,5);
    window_max=0;
    last_qrs=0;
    
    %% Filter Stage
    %figure;    
    
    % Lowpass Filter 
    ecg_2=filter(b_low,1,ecg);
    %plot(ecg_2);
 
    % Highpass Filter 
    ecg_3=filter(b_high,1,ecg_2);

    %plot(ecg_3);

    
    % Subtract mean
    ecg_4=ecg_3-mean(ecg_3);
    %plot(ecg_4)

    
    % Absolute
    ecg_5=abs(ecg_4);
    %plot(ecg_5)

    
    % Average
    ecg_6=filter(b_avg,1,ecg_5);
    %plot(ecg_6);

    
    %% Detection
    qrs_loc=zeros(length(ecg),1);
    
    ecg_res=circshift(ecg_6,[0 -round(delay)]);
    
    for (i=1:length(ecg_res))

        % Check for end of window and adjust threshold
        if (mod(i,window)==0)
            %Buffer shift
            window_max_buff=circshift(window_max_buff,[0 1]);
            window_max_buff(1) = window_max;
            h_thresh=h_thresh_correct*median(window_max_buff);
            window_max=0;
        end
        % Check for new window max
        if(ecg_res(i)>window_max)
            window_max=ecg_res(i);
        end
        % Check if refractory period is over
        if(i-last_qrs>T_ref)
            % Check if a candidate QRS was detected
            if candidate_detected==1
                % if end of candidate search was reached
                if i==end_cand_search
                    qrs_loc(candidate_pos)=1;
                    last_qrs=candidate_pos;
                    candidate_detected=0;
                    candidate_pos=0;
                    candidate=0;
                % Check if new candidate is found
                % Candidate QRS value is the maximum value and position
                % from initial high threshold crossing to a refractory
                % period after that
                elseif ecg_res(i)>candidate
                    candidate=ecg_res(i);
                    candidate_pos=i; 
                end

            else
                % Check if high threshold is surpassed 
                if(ecg_res(i)>h_thresh)
                    % Make this position the first candidate value
                    candidate=ecg_res(i);
                    candidate_pos=i;
                    candidate_detected=1;
                    end_cand_search=i+T_ref;
                end            
            end

        end
    end
    qrs=qrs_loc;
end