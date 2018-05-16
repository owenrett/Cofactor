def calc_pm_new(Pm, cm):
    import numpy as np
    msz = Pm.shape[2]
    Pm_new = np.zeros(shape = (3,3,msz))

    for i in range(Pm.shape[2]):
        Pm_new[:,:,i] = np.matmul(cm, Pm[:,:,i])

    return Pm_new
