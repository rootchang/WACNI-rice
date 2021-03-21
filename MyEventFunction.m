function [VALUE, ISTERMINAL, DIRECTION] = MyEventFunction(~,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MyEventFunction.m: This function catches exceptional events during calling of simulation_main.m by RunModel.m.
%    Currently we use this function to restrict running time of simulation_main.m to be <20 min.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TimeOut = 20*60; % 20min
    VALUE = toc-TimeOut;
    ISTERMINAL=1;
    DIRECTION = 0;
    if VALUE>0
        errorStruct.message='Time out. Possibily singular values';
        errorStruct.identifier='eCM:Timeout';
        error(errorStruct);
    end
end