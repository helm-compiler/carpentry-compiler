import os
import sys
import numpy as np

time_dict = {
(0, 0):0,    
(0, 1): 0.5,
(1, 0): 0,
(0, 2): 0,
(2, 0): 0,
(0, 3): 0,
(3, 0): 0,
(1, 1):0,
(1, 2): 0.5,
(2, 1):0,
(1, 3):0,
(3, 1):0.5,
(2, 2):0,
(2, 3): 0,
(3, 2):0,
(3, 3):0,
(0, 4):0,
(4, 0):0,
(1, 4):0,
(4, 1):0,
(2, 4):0,
(4, 2):0,
(3, 4):0,
(4, 3):0,
(4, 4):0
}

scores = [(3.75, 1.02062 ,10.5), (3.85, 1.00075, 10.5), (3.75, 1.004 ,11), (3.85, 1.06287 ,10), (10.425 ,5.0593,54.5) ]
lines = [18, 17,18, 17, 10]


# ''' For three chairs '''
# output = open("output.txt",'w+');
# for p1 in range(0, 4):
#     for p2 in range(0, 4):
#         p2_fbtime = time_dict[(p1, p2)]
#         for p3 in range(0, 4):
#             fbtime = p2_fbtime + time_dict[(p2, p3)] + time_dict[(p3, 4)]
#             precision = lines[p1] * scores[p1][1] + lines[p2] * scores[p2][1] + lines[p3] * scores[p3][1] + lines[4]* scores[4][1]
#             precision = precision / (lines[p1] + lines[p2] + lines[p3] + lines[4])
#             fin = np.asarray(scores[p1]) + np.asarray(scores[p2]) + np.asarray(scores[p3]) + np.asarray(scores[4])  
#             fin[2] -= fbtime
#             fin[1] = precision
#             output.write("%f %f %f\n" % (fin[0], fin[1], fin[2]))


# output.close()

''' For six chairs '''
output = open("output.txt",'w+');
for p1 in range(0, 4):
    for p2 in range(0, 4):
        p2_fbtime = time_dict[(p1, p2)]
        for p3 in range(0, 4):
            p3_fbtime = p2_fbtime + time_dict[(p2, p3)]
            for p4 in range(0, 4):
                p4_fbtime = p3_fbtime + time_dict[(p3, p4)]
                for p5 in range(0, 4):
                    p5_fbtime = p4_fbtime + time_dict[(p4, p5)]
                    for p6 in range(0, 4):
                        fbtime = p5_fbtime + time_dict[(p5, p6)] + time_dict[(p6, 4)]
                        precision = lines[p1] * scores[p1][1] + lines[p2] * scores[p2][1] + lines[p3] * scores[p3][1] + lines[p4] * scores[p4][1] + lines[p5]*scores[p5][1] + lines[p6]*scores[p6][1] +lines[4]* scores[4][1]
                        precision = precision / (lines[p1] + lines[p2] + lines[p3] + lines[p4] + lines[p5] +lines[p6] + lines[4])
                        fin = np.asarray(scores[p1]) + np.asarray(scores[p2]) + + np.asarray(scores[p3]) + np.asarray(scores[p4])+ np.asarray(scores[p5])+ np.asarray(scores[p6]) + np.asarray(scores[4])   
                        fin[2] -= fbtime
                        fin[1] = precision
                        output.write("%f %f %f\n" % (fin[0], fin[1], fin[2]))


output.close()