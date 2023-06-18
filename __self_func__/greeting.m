function [SubID,SubName,SubAge,SubGender,Handedness] = greeting()
Prompt = {'Subject Number', 'Name', 'Age', 'Gender:(1 is for male, 2 for female)', 'Handedness(L/R)'};
DlgTitle = 'Personalia';
Numlines = 1;
Defaultans = {'00', 'john','20','1','R'};
Answer = inputdlg(Prompt, DlgTitle, Numlines, Defaultans);
SubID = str2double(Answer{1});
SubName = Answer{2};
SubAge = str2double(Answer{3});
SubGender = str2double(Answer{4});
Handedness = Answer{5};
end
