clc;
clear;
%% random && initial
r = randperm(65); 
r = sort(r(1:35));
r(36) = 0;

test_set = {};
test_folder = [];
test_integer = 1;
train_set = {};
train_folder = [];
train_integer = 1;

%% read all image & separate into test and train
for i = 1:39
    if i == 14  % no 14th person
        continue;
    end
    dirs = pwd();
    path = fullfile(dirs, 'CroppedYale', sprintf('yaleB%02d', i));
    img_path_list = [];
    img_path_list = dir(fullfile(path, '*.pgm'));
    img_path_list = [img_path_list ; dir(fullfile(path, '*.bad'))];
    r_index = 1;
    for j = 1:65
        if strfind(img_path_list(j).name, 'Ambient') % Ambient --> 480 * 640 delete
            continue
        end
        if strfind(img_path_list(j).name, 'bad')
            continue
        end
        imdata = imread(fullfile(path, img_path_list(j).name)); % 192*168 unit8
        if j == r(r_index)
           train_set{train_integer} = double(imdata); 
           train_folder = [train_folder, i];
           r_index = r_index + 1;
           train_integer = train_integer + 1;
           continue
        end
        test_set{test_integer} = double(imdata); 
        test_folder = [test_folder, i];
        test_integer = test_integer + 1;
    end
end

%% use "SSD"(Sum of Squared Differences), "SAD"(Sum of Absolute Difference) to find & calculate accuracy
[row, test_count] = size(test_set);
[row, train_count] = size(train_set);
SAD_folder = -1;
SAD_correct = 0;
SSD_folder = -1;
SSD_correct = 0;
for i = 1:test_count
    SAD_min = Inf;
    SSD_min = Inf;
    for j = 1:train_count  
        value = test_set{i} - train_set{j};
        %SAD
        diff_SAD_matrix =sum(abs(value), "all") ;
        if diff_SAD_matrix < SAD_min
            SAD_min = diff_SAD_matrix;
            SAD_folder = train_folder(j);
        end
        %SSD
        diff_SSD_matrix = sum(value(:).^2);
        if diff_SSD_matrix < SSD_min
            SSD_min = diff_SSD_matrix;
            SSD_folder = train_folder(j);
        end
    end
    if test_folder(i) == SAD_folder
        SAD_correct = SAD_correct + 1;
    end
    if test_folder(i) == SSD_folder
        SSD_correct = SSD_correct + 1;
    end
end

%% printf
fprintf('SAD : %f%%\nSSD : %f%%\n', SAD_correct/test_count*100, SSD_correct/test_count*100);
