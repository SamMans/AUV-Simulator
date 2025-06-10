         % Author: Samer A. Mohamed (https://github.com/SamMans) %
classdef Thruster
    properties(SetAccess = private)
        ID %--> Thruster ID
        B_table %--> Bluerobotics table of thrust
        F %--> Current thrust force in Newtons
    end
    methods
        function obj = Thruster(i, model, volt)
            % Identify thruster number
            obj.ID = i; 

            % Read the table of thrust from excel sheet
            obj.B_table = readmatrix(pwd + "\Libraries\Thrusters\" + model + ".xls", 'Sheet', volt);
        end
        function obj = P2T(obj, pwm) %--> PWM to Thrust converter
            % Thrust index/location in table
            [~,I] = min(abs(pwm - obj.B_table(:,1)));
            
            % Return actual thrust + Noise
            obj.F = obj.B_table(I, 6) * 9.80665; 
        end
        function tt = getT(obj)
            % Return the current thrust value
            tt = obj.F; 
        end
    end
end