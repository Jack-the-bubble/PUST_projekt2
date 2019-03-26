function MinimalWorkingExample()
    addpath('F:\SerialCommunication'); % add a path to the functions
    initSerialControl COM5 % initialise com port

    step_response=[];
    figure;
    while(1)
        %% obtaining measurements
        measurements = readMeasurements(1:7); % read measurements from 1 to 7
        
        %% processing of the measurements and new control values calculation
        disp(measurements); % process measurements
         step_response=[step_response measurements(1)];
         plot(step_response);
         %save('step-response28z60');
   
        %% sending new values of control signals
        %sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
        %             [50, 0, 0, 0, 28, 0]);  % new corresponding control values
        sendControlsToG1AndDisturbance(28,0);
        
        %% synchronising with the control process
        waitForNewIteration(); % wait for new batch of measurements to be ready
        drawnow;
    end
end