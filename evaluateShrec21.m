clear all;

results = {};

fid = fopen('results.txt');
for i=1:72
    tmpLine = fgetl(fid);
    tmpArray = strsplit(tmpLine,';');
    results{i,1} = str2num(tmpArray{1,1});
    results{i,2} = tmpArray(1,2:end-1);
end
fclose(fid);

results = sortrows(results,1);

annotations = {};

fid = fopen('annotations.txt');
for i=1:72
    tmpLine = fgetl(fid);
    tmpArray = strsplit(tmpLine,';');
    annotations{i,1} = str2num(tmpArray{1,1});
    annotations{i,2} = tmpArray(1,2:end-1);
end
fclose(fid);

annotations = sortrows(annotations,1);

numGest = 18;
gestures = ["ONE" "TWO" "THREE" "FOUR" "OK" "MENU" "POINTING" "LEFT" "RIGHT" "CIRCLE" "V" "CROSS" "GRAB" "PINCH" "TAP" "DENY" "KNOB" "EXPAND"];

minOverlapRatio = 0.5;
%classResults columns: Total, Correct, Missed, Misclassified, FalsePositive, Overlap
classResults = zeros(numGest,6);
classResults(:,1) = zeros(numGest,1)+16;

for s = 1:72
A = annotations{s,2};
R = results{s,2};
found = zeros(1,size(A,2)/3);
for r = 1:3:size(R,2)
    RR = R(1,r:r+2);
    detected = false;
    countA = 0;
    for a = 1:3:size(A,2)
        countA=countA+1;
        AA = A(1,a:a+2);
        
        AAlength = str2num(AA{3})-str2num(AA{2});
        %RRlength = str2num(RR{3})-str2num(RR{2});
        overlap = min([str2num(AA{3}) str2num(RR{3})])-max([str2num(AA{2}) str2num(RR{2})]);
        overlapRatio = overlap/AAlength;
        
        if (overlapRatio>minOverlapRatio) %something has been detected where a labeled gesture exists
            detected = true;
            if (strcmp(RR{1},AA{1})) % the lablel matches, it's a correct detection
                if (found(1,countA)~=1) %to avoid counting more than once a duplicated detection 
                    found(1,countA) = 1; %marking that a correct dection for this gesture has been found
                    classResults(find(gestures==AA{1}),2) = classResults(find(gestures==AA{1}),2)+1;
                    classResults(find(gestures==AA{1}),6) = classResults(find(gestures==AA{1}),6)+overlapRatio;
                end
            else %the label doesn't  match, it's a misclassification
                classResults(find(gestures==AA{1}),4) = classResults(find(gestures==AA{1}),4)+1;
            end         
        end
        if(detected)
            break;
        end
        
    end
    if (~detected) %the detected gesture window did not match with any of the labeled one, it's a false positive
        classResults(find(gestures==AA{1}),5) = classResults(find(gestures==AA{1}),5)+1;
    end
    
end
for f = 1:size(found,2)
    if (found(1,f)==0) %this labeled gesture was neither detected correctly or misclassified
        classResults(find(gestures==A{1+((f-1)*3)}),2) = classResults(find(gestures==A{1+((f-1)*3)}),2)+1;
    end
end
end
classResults(:,6) = classResults(:,6)./classResults(:,2); %averaging the overlapping ratio of the gesture windows (detected vs. labeled)

classPrecision = classResults(:,2)./(classResults(:,2)+classResults(:,3)+classResults(:,5));
classRecall = classResults(:,2)./(classResults(:,2)+classResults(:,4));

correctScore = sum(classResults(:,2))/sum(classResults(:,1))
misclassifiedRate = sum(classResults(:,4))/sum(classResults(:,1))
falspositiveRate = sum(classResults(:,5))/sum(classResults(:,1))
precision = sum(classPrecision)/numGest
recall = sum(classRecall)/numGest















