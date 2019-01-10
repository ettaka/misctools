import numpy as np
import itertools
import cpickle as pickle

def create_harmonics_data():
    try:
        print "tot_harmonics.pickle detected, using that"
        tot_harmonics = pickle.load("tot_harmonics.pickle")
    except:
        print "tot_harmonics.pickle not detected, creating data"
        s_array = []
        s_array += [np.array([[ 1,-1],
                        [ 1,-1],
                        [ 1,-1]])]
        s_array += [np.array([[ 1, 1],
                        [ 1, 1],
                        [ 1, 1]])]
        s_array += [np.array([[ 1, 1],
                        [-1, 1],
                        [-1,-1]])]
        s_array += [np.array([[-1, 1],
                        [-1,-1],
                        [ 1,-1]])]
        s_array += [np.array([[-1, 1],
                        [ 1,-1],
                        [-1, 1]])]
        s_array += [np.array([[-1,-1],
                        [ 1, 1],
                        [-1,-1]])]
        s_array += [np.array([[-1,-1],
                        [-1, 1],
                        [ 1, 1]])]
        s_array += [np.array([[ 1,-1],
                        [-1,-1],
                        [-1, 1]])]

        s_array = np.array(s_array)

        s1_harmonics = np.array([[ 2.1,-0.49],
                                 [0.73, 0.21],
                                 [0.19, 0.10]])

        MQXFS6a_harmonics = np.array([[-4.26,  3.06],
                                      [ 0.66, -1.59],
                                      [-0.34, -1.39]])

        harmonics = MQXFS6a_harmonics

        combination_vector = (0,1,2,3,4,5,6,7)
        all_combinations = []
        for comb_element in combination_vector:
            all_combinations += list(itertools.combinations(combination_vector,comb_element+1))

        tot_harmonics = []
        for combination in all_combinations:
            comb = np.array(combination)
            delta_harmonics = s1_harmonics * np.sum(s_array[comb],0)
            tot_harmonics.append(harmonics + delta_harmonics)

        print "Dump tot_harmonics.pickle"
        pickle.dump(tot_harmonics, "tot_harmonics.pickle")

    return tot_harmonics

tot_harmonics = create_harmonics_data()

print np.abs(tot_harmonics)[:,0,:]

    


