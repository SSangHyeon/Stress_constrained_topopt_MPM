function importfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  지정된 파일에서 데이터 가져오기
%  FILETOREAD1:  읽을 파일

%  MATLAB에서 07-Feb-2022 10:00:04에 자동 생성됨

% 파일 가져오기
newData1 = load('-mat', fileToRead1);

% 기본 작업 공간에 해당 필드로부터 새 변수를 생성합니다.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('caller', vars{i}, newData1.(vars{i}));
end

