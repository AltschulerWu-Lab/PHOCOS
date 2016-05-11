function [Precision,Recall,F]=CalculatePecision(TP,FP,TN,FN)

Precision=TP/(TP+FP);
Recall=TP/(TP+FN);

F=2/(1/Precision+1/Recall);