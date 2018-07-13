function score = o_transferScore(PowerdB,mode)
% TRANSFERSCORE - Computes a score that reflects how good the
% directivity value is in the MU-MIMO scenario. Directivities start
% getting an score when their absolute value is greater than 0.
% Directivities are capped up at the maximum tolerable antenna gain by
% the FCC community (33dB if Ptx=10dBm)
%
% For more information, see slide 7-4 in:
% https://www.cse.wustl.edu/~jain/cse574-14/ftp/j_07sgh.pdf
%
% Syntax:  score = transferScore(PowerdB,mode)
%
% Inputs:
%    PowerdB - Description
%    mode - 1 for Directivity to intended user. 0 to compute
%    directivity towards interfeered users.
%
% Outputs:
%    score - Value between 0 and 1 that reflects how good the current
%    antenna selection is.
%
% Example: 
%    score = transferScore(10,1)  % Compute score for intended user
%    score = transferScore(-10,0)  % Compute score to other users
%
%------------- BEGIN CODE --------------

% Previously used
%     % Define the scoring function here. Left hand-side values reflect the
%     % power levels in dB of the directivity. The right hand-side values
%     % reflect the score obtained
%     p = [ -500 -1       1;
%            -50 -1       1;
%            -40 -1     0.9;
%            -33 -1    0.85;
%            -30 -0.9   0.8;
%            -25 -0.8   0.7;
%            -20 -0.6   0.6;
%            -15 -0.4   0.4;
%            -10 -0.3   0.3;
%             -5 -0.2   0.2;
%              0 0        0;
%              5 0.2   -0.2;
%             10 0.3   -0.3;
%             15 0.4   -0.4;
%             20 0.6   -0.6;
%             25 0.8   -0.7;
%             30 0.9   -0.8;
%             33 1    -0.85;
%             40 1     -0.9;
%             50 1       -1;
%            500 1       -1].';
% 	x = p(1,:);  % Power in dB
%     if mode
%         % Prx (intended user)
%         y = p(2,:);  % Score
%     else
%         % Int (interfereed users)
%         y = p(3,:);  % Score
%     end
%     xx = min(x):0.05:max(x);
%     yy = interpn(x,y,xx,'linear');
%     [~,idx] = min(abs(PowerdB - xx));
%     score = yy(idx);

% Alternative (linear)
 if mode    
     mindB = -33;
     maxdB = 33;
     minScore = -1;
     maxScore = 1;
     m = (maxScore-minScore)/(maxdB-mindB);  % Slope 
     x = (mindB:0.01:maxdB);
     n = maxScore - m*maxdB;
     y = x.*m + n;
     % Include limits
     x = [-500 x 500];
     y = [-1 y 1];
 else
     % 1st stage
     mindB1st = 50;
     maxdB1st = -50;
     minScore1st = -1;
     maxScore1st = 0.85;
     m = (maxScore1st-minScore1st)/(maxdB1st-mindB1st);  % Slope 
     x1 = (maxdB1st:0.01:mindB1st);
     n1 = maxScore1st - m*maxdB1st;
     y1 = x1.*m + n1;
     % 2nd stage
     mindB2nd = maxdB1st;
     maxdB2nd = -100;
     minScore2nd = maxScore1st;
     maxScore2nd = 1;
     m = (maxScore2nd-minScore2nd)/(maxdB2nd-mindB2nd);  % Slope 
     x2 = (maxdB2nd:0.01:mindB2nd);
     n2 = maxScore2nd - m*maxdB2nd;
     y2 = x2.*m + n2;
     % Append stages
     x = [x2 x1];
     y = [y2 y1];
     % Include limits
     x = [-500 x 500];
     y = [1 y -1];
 end
 [~,idx] = min(abs(PowerdB - x));
 score = y(idx);

     
     
% EOF