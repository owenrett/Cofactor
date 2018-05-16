def axis_gener(A, B):
    import numpy as np
    #generates rotation axis, e, from input matrices A and B

    #assuming QB - A = a \items n

    C = np.matmul(np.linalg.inv(A), np.matmul(np.matmul(B, B), np.linalg.inv(A)))

    #sorting algorithm from https://stackoverflow.com/questions/8092920/sort-eigenvalues-and-associated-eigenvectors-after-using-numpy-linalg-eig-in-pyt
    eigenValues, eigenVectors = np.linalg.eig(C)

    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]

    e1 = eigenVectors[:,2]

    e2 = eigenVectors[:, 1]

    e3 = eigenVectors[:, 0]

    la1 = eigenValues[2]

    la2 = eigenValues[1]

    la3 = eigenValues[0]

    s = - (np.dot(e2, np.matmul(np.linalg.matrix_power(A,2) , e1)) / (np.sqrt(la3) * np.dot(e2, np.matmul(np.linalg.matrix_power(A,2) , e3))))

    if s < 0:
        s = -1
    else:
        s = 1

    del1 = np.dot(e1, np.matmul(np.linalg.matrix_power(A,2) , e1)) + s * np.sqrt(la3) * np.dot(e3, np.matmul(np.linalg.matrix_power(A,2) , e1))

    del1 = 1 / np.sqrt(2 * np.abs(del1))

    del3 = s * np.sqrt(la3) * del1

    e = (del1 * np.matmul(A ,e1) + del3 * np.matmul(A , e3))

    e = e / np.linalg.norm(e)

    e[np.abs(e) < 10^-4] = 0

    if np.isnan(e).any():
        e = np.array([[0], [0], [0]])

    return e