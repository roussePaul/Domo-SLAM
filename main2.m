path(path,'../test_cases_lab1_v0.2');

test_cases_lab1();

%%
path(path,'../Lab2/Datasets');

runlocalization_track('so_o3_ie.txt', 'map_o3.txt', 1, 1, 1, 2);


%%

runlocalization_track('so_pb_10_outlier.txt', 'map_pent_big_10.txt', 1, 1, 1, 2);

%%

runlocalization_track('so_pb_40_no.txt', 'map_pent_big_40.txt', 1, 1, 1, 2);
