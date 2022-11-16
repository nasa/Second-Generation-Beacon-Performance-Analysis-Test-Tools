%
%
%   function str2bin(string)
%
%   converts a string of binary digits into a vector of binary bits
%
%
function o=str2bin(x)

for i=1:length(x)
    if(x(i) == '1')
        o(i) = 1;
    else
        o(i) = 0;
    end
end