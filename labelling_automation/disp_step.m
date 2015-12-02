function disp_step(msg)
% displays the message with the step number and a newline
% eg: 1) Do something \n
% 
% the step increments by 1 each time the method is called
% 
% call the method with an empty string to print nothing and reset the step
% counter, eg disp_step('');

    persistent step;
    
    if (isempty(step))
        step = 1;
    else
        step = step + 1;
    end
    
    if (isempty(msg))
        step = 0;
        return;
    end
    
    fprintf(['%d) ' msg '\n'], step);
end