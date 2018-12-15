import sys

idxs = [(1, 2, 1),
 (1, 4, 2),
 (2, 2, 1),
 (4, 4, 1),
 (5, 1, 3),
 (7, 3, 1),
 (8, 4, 1),
 (21, 2, 2),
 (21, 4, 1),
 (22, 4, 1),
 (24, 1, 1),
 (25, 1, 2),
 (25, 4, 1),
 (27, 2, 1),
 (29, 4, 1),
 (30, 4, 1),
 (31, 2, 1),
 (31, 3, 2),
 (32, 2, 1),
 (32, 4, 1),
 (33, 3, 2),
 (34, 2, 1),
 (35, 2, 1),
 (36, 2, 1),
 (38, 2, 1),
 (39, 4, 1),
 (40, 2, 1),
 (41, 1, 2),
 (41, 2, 1),
 (41, 3, 1),
 (42, 2, 4),
 (45, 1, 1),
 (48, 4, 1),
 (54, 2, 1),
 (54, 3, 1),
 (54, 4, 1),
 (56, 3, 1),
 (58, 2, 1),
 (59, 2, 2),
 (59, 4, 1),
 (60, 4, 1),
 (61, 1, 1),
 (61, 4, 2),
 (62, 3, 1),
 (65, 4, 1),
 (66, 1, 1),
 (71, 1, 1),
 (71, 3, 1),
 (72, 4, 1),
 (74, 1, 1),
 (75, 4, 1),
 (76, 1, 3),
 (81, 4, 1),
 (82, 2, 1),
 (82, 4, 1),
 (83, 3, 1),
 (85, 4, 2),
 (87, 1, 2),
 (87, 2, 1),
 (88, 1, 1),
 (88, 4, 1),
 (91, 2, 1),
 (94, 2, 1),
 (97, 4, 1),
 (98, 2, 2),
 (101, 1, 1),
 (101, 2, 1),
 (102, 1, 1),
 (102, 4, 1),
 (103, 1, 1),
 (103, 4, 1),
 (104, 2, 1),
 (104, 4, 1),
 (105, 4, 2),
 (108, 2, 2),
 (108, 4, 1),
 (111, 3, 1),
 (112, 2, 1),
 (112, 4, 1),
 (116, 1, 1),
 (116, 3, 1),
 (116, 4, 2),
 (117, 1, 1),
 (118, 4, 1),
 (119, 2, 1),
 (122, 2, 3),
 (123, 1, 1),
 (125, 2, 1),
 (126, 1, 1),
 (126, 3, 1),
 (128, 2, 2),
 (131, 1, 1),
 (131, 2, 1),
 (131, 4, 1),
 (132, 4, 1),
 (135, 2, 1),
 (136, 1, 1),
 (138, 4, 1),
 (141, 1, 1),
 (141, 3, 1),
 (147, 4, 1),
 (148, 2, 1),
 (148, 4, 1),
 (149, 1, 1),
 (149, 4, 1),
 (150, 2, 1),
 (151, 1, 1),
 (152, 2, 1),
 (153, 2, 1),
 (153, 4, 1),
 (156, 4, 1),
 (157, 2, 1),
 (158, 1, 1),
 (160, 2, 1),
 (161, 1, 1),
 (161, 2, 2),
 (161, 3, 2),
 (162, 1, 1),
 (164, 1, 1),
 (165, 1, 1),
 (165, 4, 1),
 (167, 1, 1),
 (167, 2, 1),
 (169, 2, 1),
 (169, 4, 1),
 (171, 1, 2),
 (171, 3, 1),
 (175, 2, 1),
 (175, 4, 1),
 (176, 2, 1),
 (178, 1, 2),
 (180, 4, 2),
 (181, 2, 1),
 (182, 4, 2),
 (183, 1, 1),
 (183, 3, 1),
 (184, 2, 1),
 (187, 4, 1),
 (188, 4, 1),
 (191, 1, 2),
 (192, 1, 1),
 (192, 4, 1),
 (193, 4, 1),
 (194, 4, 1),
 (197, 2, 1)]

# Get datafile and option
# opt = 0 is the dfile
# opt = 1 is the model
# opt = 2 is the number of iterations
idx = int(sys.argv[1])
opt = int(sys.argv[2])

# Print requested result for bash output
print(idxs[idx][opt])
