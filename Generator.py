def generator(cell_name):
    import numpy as np
    #this is messy beyond hell, but just returns all symmetry operation matrices related to the given cell name
    #calling UC interp to get generators for the unit cell
    #programming in a very limited ability to interpret non-string names, but better safe than sorry
    #using string names of cell, ex. generator("m-3m")

    holding = UC_interp(cell_name)

    k = 0

    a = holding.shape

    if len(a) < 3:
        a = 1
    else:
        a = holding.shape[2]

    while k <= a-1:
        m = 0
        while m <= a-1:
        #iterating across all possible products of generator matrices
            if a == 1:
                temp = np.matmul(holding[:,:],holding[:,:])
            else:
                temp = np.matmul(holding[:, :, k], holding[:, :, m])
            #checking for uniqueness
            su = 0
            for n in range(a):
                if a == 1:
                    if np.array_equal(temp, holding[:, :]) or np.array_equal(temp, np.zeros(shape=(3, 3))):
                        su += 1
                else:
                    if np.array_equal(temp, holding[:, :, n]) or np.array_equal(temp, np.zeros(shape=(3, 3))):
                        su += 1

            if su == 0:
                holding = np.dstack((holding, temp))
            m = m + 1

            a = holding.shape

            if len(a) < 3:
                a = 1
            else:
                a = holding.shape[2]
        k = k + 1


    return(holding)
