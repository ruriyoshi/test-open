function [F] = method_pseudo(W,P)
Z = pinv(W);%P��F�ϊ��s��
F = Z*P;
F = method_goto0(F);
end