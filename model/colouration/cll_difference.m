function [CLL_difference,freqs] = cll_difference(input,ref,Fs)
if length(input(:,1)) < length(input(1,:)) 
    input = input';
end
if length(ref(:,1)) < length(ref(1,:)) 
    ref = ref';
end

[CLL_input,freqs] = cll(input(:,1),input(:,2),Fs); 
[CLL_ref,~] = cll(ref(:,1),ref(:,2),Fs);

CLL_difference = CLL_input-CLL_ref;
end

