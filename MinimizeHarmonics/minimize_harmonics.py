import numpy as np
import itertools
import cPickle as pickle
ref_harmonics = []
#ref_harmonics += [np.array([[ 2.1 , 0.49],
#                      [ 0.73,-0.21],
#                      [ 0.19,-0.1]])]
#
#ref_harmonics += [np.array([[ 2.1 ,-0.49],
#                      [ 0.73, 0.21],
#                      [ 0.19, 0.1]])]
#
#ref_harmonics += [np.array([[-0.49, 2.1],
#                      [-0.73, 0.21],
#                      [-0.1 ,-0.19]])]
#
#ref_harmonics += [np.array([[ 0.49, 2.1],
#                      [-0.73,-0.21],
#                      [ 0.1 ,-0.19]])]
#
#ref_harmonics += [np.array([[-2.1 ,-0.49],
#                      [ 0.73,-0.21],
#                      [-0.19, 0.1]])]
#
#ref_harmonics += [np.array([[-2.1 , 0.49],
#                      [ 0.73, 0.21],
#                      [-0.19,-0.1]])]
#
#ref_harmonics += [np.array([[ 0.49,-2.1],
#                      [-0.73, 0.21],
#                      [ 0.1 , 0.19]])]
#
#ref_harmonics += [np.array([[-0.49,-2.1],
#                      [-0.73,-0.21],
#                      [-0.1 , 0.19]])]

ref_harmonics += [np.array([[ 1.80, 0.11 ],
                            [ 0.59,-0.26 ],
                            [ 0.17,-0.11 ]])]

ref_harmonics += [np.array([[ 1.80,-0.11 ],
                            [ 0.59, 0.26 ],
                            [ 0.17, 0.11 ]])]

ref_harmonics += [np.array([[-0.11, 1.80 ],
                            [-0.59, 0.26 ],
                            [-0.11,-0.17 ]])]

ref_harmonics += [np.array([[ 0.11, 1.80 ],
                            [-0.59,-0.26 ],
                            [ 0.11,-0.17 ]])]

ref_harmonics += [np.array([[-1.80,-0.11 ],
                            [ 0.59,-0.26 ],
                            [-0.17, 0.11 ]])]

ref_harmonics += [np.array([[-1.80, 0.11 ],
                            [ 0.59, 0.26 ],
                            [-0.17,-0.11 ]])]

ref_harmonics += [np.array([[ 0.11,-1.80 ],
                            [-0.59, 0.26 ],
                            [ 0.11, 0.17 ]])]

ref_harmonics += [np.array([[-0.11,-1.80 ],
                            [-0.59,-0.26 ],
                            [-0.11, 0.17 ]])]

ref_harmonics = np.array(ref_harmonics)

s2_harmonics = np.array([[ 2.1,-0.49],
                         [0.73, 0.21],
                         [0.19, 0.10]])

MQXFS6a_harmonics = np.array([[-4.26,  3.06],
                              [ 0.66, -1.59],
                              [-1.39, -1.39]])

def compute_delta_harmonics(ref_harmonics, combination):
    comb=list(combination)
    delta_harmonics =  np.sum(ref_harmonics[comb],0)
    #delta_harmonics =  np.sum(ref_harmonics[comb],0)
    #delta_harmonics =  ref_harmonics[comb]
    return delta_harmonics

def create_harmonics_data():
    harmonics = MQXFS6a_harmonics

    combination_vector = (0,1,2,3,4,5,6,7)
    all_combinations = []
    for comb_element in combination_vector:
        all_combinations += list(itertools.combinations(combination_vector,comb_element+1))

    tot_harmonics = []
    for combination in all_combinations:
        comb = np.array(combination)
        delta_harmonics = compute_delta_harmonics(ref_harmonics, combination)
        tot_harmonics.append(harmonics + delta_harmonics)
    tot_harmonics = np.array(tot_harmonics)

    return tot_harmonics, all_combinations 

#print compute_delta_harmonics(ref_harmonics, (0,1,6,7))
#exit()

number_of_printed_combinations = 10
min_row = 0
tot_harmonics, all_combinations = create_harmonics_data()
#min_row_sums = np.sum(np.abs(tot_harmonics)[:,min_row,:],1) # minimum row sum absolute
#min_row_sums = np.abs(np.sum(tot_harmonics[:,min_row,:],1)) # minimum row absolute sum
min_row_sums = np.sum(np.sum(np.abs(tot_harmonics)[:,:,:],1),1)# min sum sum absolute
sorted_list_indices = np.argsort(min_row_sums)
for i in range(number_of_printed_combinations):
    combination_index = sorted_list_indices[i]
    print "__________________________"
    print "combination_index:", combination_index, "combination:", np.array(all_combinations[combination_index])+1
    print "total harmonics:\n", tot_harmonics[combination_index]#.flatten()



    


