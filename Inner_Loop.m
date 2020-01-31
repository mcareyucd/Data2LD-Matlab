function [ output_args ] = Inner_Loop( in1, in2, in3, in4 ,wvec )
%InnerLoop Summary of this function goes here
%   Detailed explanation goes here

output_args = inner_loop(in1', in2', in3', in4', wvec)';
% output_args = loop_mex(in1', in2', in3', in4');

end

