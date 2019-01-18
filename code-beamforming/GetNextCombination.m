%[A,M,H,ITERATION] = GETNEXTCOMBINATION(N, K, A, M, H, ITERATION) returns 
%the next K-subset of an N-set. It implements the "Revolving Door" (NEXKSB)
%algorithm found in "Combinatorial Algorithms", 2nd Ed., by Albert Nijenhuis
%and Herbert S. Wilf. It is a recursive alternative to the NCHOOSEK command
%that is extremely fast (O(1)), at the expense of slightly more user
%interaction. (You do get something that NCHOOSEK does not give you for this
%extra effort: see input "m", below.) Apart from its efficiency, it is
%awesome because each new K-subset differs from the last by a single element
%and the last combination generated differs from the first combination by a
%single element. 
%
% Inputs: N - a natural number corresponding to the length of the N-set
%         K - a natural number, less than N, corresponding to the length of
%             the K-subset
%         A - the current K-subset                 
%         m - index used internally to generate next combination. It must be 
%             initialized to 0 by the user and that is all the interaction 
%             needed.
%             ***m can be used to track which element came out of A and
%             which came in. E.g., for the current combination A, m came out
%             of the previous combination and was replaced by (m+1).
%         h - index used internally to generate next combination. It must be 
%             initialized to K by the user and that is all the interaction 
%             needed.
% iteration - a natural number that counts the number of combinations. When 
%             the last pattern is generated, it will become negative.
%
% WARNINGS: 
% 1. A(N,K), with K=0, is defined as empty. It doesn't make sense to add code 
%    to handle this one case; so we demand that K be a natural number, i.e. K>=1.
% 2. The initial A is assumed to be ordered as follows:
%    1<=A(1)<A(2)<...<=N. This doesn't dimish the utility of the algorithm
%    at all. You can always set A = 1:K.
% 
% Example: (Replaces "combinations = nchoosek(V,K)")
%
%         V = [4,1,3,5,7];  %THIS IS YOUR VECTOR
%         N = 5;            %LENGTH OF V
%         K = 3;           
%         m = 0;
%         h = K;
%         iteration = 1;         
%         A = 1:K;          %FIRST COMBINATION
%         patterns = [];    %YOU CAN COLLECT THEM, BUT IT'S SLOWER
%         while(iteration>0)
%             [A,m,h,iteration] = GetNextCombination(N, K, A, m, h, iteration);
%             patterns = [patterns;V(A)];
%         end
%
% COMMENTS: GETNEXTCOMBINATION really shines when you DO NOT collect the
%           combinations as in the example above, and when nchoosek(N,K) is
%           a huge number. 
%
%           I've left comments in the code in case you want to compare it
%           to the original FORTRAN code. I've replaced IFs with
%           multiplications because it's faster and FOR loops with colon
%           indexing.
%       
% REFERENCE: "Combinatorial Algorithms", 2nd Ed., by Albert Nijenhuis and 
%             Herbert S. Wilf.
%             http://www.math.upenn.edu/~wilf/website/CombAlgDownld.html
%
% Script written by: Jose A. Lopez
%                    IGERT Fellow
%                    Robust Systems Lab
%                    Northeastern University
%
%This code is licensed under the BSD License as required by MATLAB Central.
%
%Live long and prosper.
%
function [A,m,h,iteration,out] = GetNextCombination(N, K, A, m, h, iteration)
if(iteration~=1)
%     if(m<(N-h))
%         h = 0;
%     end
    h = (m>=(N-h))*h;
    h = h + 1;    
    m = A(K+1-h);
else
    m = 0;
    h = K;
end

A(K+(1:h)-h) = m + (1:h);
% if(A(1)~=(N-K+1))
%     iteration = iteration + 1;
% else
%     iteration = -iteration;
% end
iteration = (A(1)~=(N-K+1))*(iteration + 1)-(A(1)==(N-K+1))*iteration;
