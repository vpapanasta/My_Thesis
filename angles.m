% This script calculates  mean and standard deviation for the angles being 
% saved in meth_angles.mat and prop_angles.mat files

clc;
meth = load('meth_angles.mat');

disp('*****************************');
disp('*****************************');
for f = fieldnames(meth)'
   disp(['Field named: ' f{1} ]);
   disp('S.T.D: ')
   disp(std(meth.(f{1})));
   disp('MEAN: ')
   disp(mean(meth.(f{1})));
end

clc;
prop = load('prop_angles.mat');

disp('*****************************');
disp('*****************************');
for f = fieldnames(prop)'
   disp(['Name: ' f{1} ]);
   disp('S.T.D: ')
   disp(std(prop.(f{1})));
   disp('MEAN: ')
   disp(mean(prop.(f{1})));
end