%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������O��web���烍�[�J���ɕۑ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = save_log()
%�������O��web����xlsx�`���Ń_�E�����[�h
ID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
url_name = strcat('https://docs.google.com/spreadsheets/d/',ID,'/export?format=xlsx');
websave('exp_log.xlsx',url_name);

