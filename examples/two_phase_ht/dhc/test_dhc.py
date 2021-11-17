def check():
    import numpy as np
    dataOK = np.loadtxt('nusselt_ref.out')
    dataChk= np.loadtxt('data/post/wall/nusselt.out')
    tol    = 1e-6
    nts    = 10000
    chk = (np.mean(dataOK[-nts:,2])-np.mean(dataChk[-nts:,2]))<tol
    return chk


def test_answer():
    assert check()
