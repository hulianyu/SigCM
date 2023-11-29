function result = ACC(Y, predY)
%if pred_classnum
res = bestMap(Y, predY);
% accuarcy
result = length(find(Y == res))/length(Y);
end