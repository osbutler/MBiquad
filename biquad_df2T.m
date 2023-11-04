classdef  biquad_df2T < handle
% biquad_Stereo_df2T - Calculate 2n order biquad coefficient
% filter (lowpass, highpass, bandpass, peak, lowshelf, highshelf). Process
% data using the direct transform 2 transposed structure.
% --------------------------
% Author:  Oscar Butler
% Project: MBiquad
% Date:    11.4.2023
% --------------------------
    
%% Properties 
    properties (Access = public)
        numStages   % Biquad order
        coeffs      % Biquad coefficient buffer
        state       % Biquad buffer state containing previous output and incoming input samples
        f0          % Significant frequency
        Q           % Quality factor of the filter
        fs          % Sampling frequency
        gaindB      % Gain in dB 
        type        % Filter type
        outputBuffer% Buffer where filtered samples are stored
    end
    methods 
        function obj = init(obj,numStages, coeffs, state, freqCut, Q, fs, gaindB, outputBuffer, type)
          obj.numStages = numStages;
          obj.coeffs = coeffs;
          obj.state = state;
          obj.f0 = freqCut;
          obj.Q = Q;
          obj.fs = fs;
          obj.gaindB = gaindB;
          obj.outputBuffer = outputBuffer;
          obj.type = type;
        end
        
        function obj = mono_df2T(obj, inputBuffer, blockSize)
            % Load coefficients
            b0 = obj.coeffs(1);
            b1 = obj.coeffs(2);
            b2 = obj.coeffs(3);
            a1 = obj.coeffs(4);
            a2 = obj.coeffs(5);

            % Load state variables
            d1L = obj.state(1);
            d2L = obj.state(2);
            d1R = obj.state(3);
            d2R = obj.state(4);

            stage = obj.numStages;

            while stage > 0
                for i=1:blockSize
                inSample = inputBuffer(i);

                outSample = b0 * inSample + d1R;
                d1R = b1 * inSample - a1 * outSample + d2R;
                d2R = b2 * inSample - a2 * outSample;

                obj.outputBuffer(i) = outSample;

                end

                % Save new state variable by incrementing state pointer
                obj.state(1) = d1L;
                obj.state(2) = d2L;
                obj.state(3) = d1R;
                obj.state(4) = d2R;

                % The current stage input is given as the output to the next stage
                inputBuffer = obj.outputBuffer;

                stage = stage - 1;
            end
        end
        
        function obj = biquad_coeff_calculation(obj)
            w0 = 2*pi*obj.f0/obj.fs;
            alpha = sin(w0)/(2*obj.Q);
            switch (obj.type)
                case 0 % LP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*(1.0 - cos(w0))/2.0;
                    b1 = A*(1.0 - cos(w0));
                    b2 = A*(1.0 - cos(w0))/2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 1 % HP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*(1.0 + cos(w0))/2.0;
                    b1 = -A*(1.0 + cos(w0));
                    b2 = A*(1.0 + cos(w0))/2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 2 % BP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*alpha;
                    b1 = 0.0;
                    b2 = -A*alpha;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 3 % PK
                    A = 10^(obj.gaindB/40.0);
                    b0 = 1.0 + alpha*A;
                    b1 = -2.0*cos(w0);
                    b2 = 1.0 - alpha*A;
                    a0 = 1.0 + alpha/A;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha/A;
                case 4 % LS
                    A = 10^(obj.gaindB/40.0);
                    b0 = A*((A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha);
                    b1 = 2.0*A*((A-1.0) - (A+1.0)*cos(w0));
                    b2 = A*((A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha);
                    a0 = (A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
                    a1 = -2.0*((A-1.0) + (A+1.0)*cos(w0));
                    a2 = (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
                case 5 % HS
                    A = 10^(obj.gaindB/40.0);
                    b0 = A*((A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha);
                    b1 = -2.0*A*((A-1.0) + (A+1.0)*cos(w0));
                    b2 = A*( (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha);
                    a0 = (A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
                    a1 = 2.0*((A-1.0) - (A+1.0)*cos(w0));
                    a2 = (A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
                otherwise
                    b0 = 1.0;
                    b1 = 0.0;
                    b2 = 0.0;
                    a0 = 1.0;
                    a1 = 0.0;
                    a2 = 0.0;
            end
            b0 = b0/a0;
            b1 = b1/a0;
            b2 = b2/a0;
            a1 = a1/a0;
            a2 = a2/a0;

            coeffs(1) = b0;
            coeffs(2) = b1;
            coeffs(3) = b2;
            coeffs(4) = a1;
            coeffs(5) = a2;
            
            obj.coeffs = coeffs;
        end
    end
end