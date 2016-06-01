function Connectivity_Ratio = connect(num_of_pre_synaptic_neurons, num_of_post_synaptic_neurons, p)
% This function will get an input of 2 numbers:
% The first being the number of pre synaptic neurons and the second being
% the number of post synaptic nerons.
% The function will genrate an mxn matrix, m being the number of
% presynaptic neurons and n being the number of post synaptic neurons.
% Also, the function will get as input the probability of connectivity. 
% Using these parameters, the function will create the matrix based on the
% binomial disturbution function.
    Connectivity_Ratio = binornd(1,p,[num_of_pre_synaptic_neurons,num_of_post_synaptic_neurons]);
end