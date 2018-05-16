def cof(matrix):
    import numpy as np
    return np.linalg.inv(matrix).T * np.linalg.det(matrix)