import numpy as np
import matplotlib.pyplot as plt

# Fig 3.3: Term structure of cash
def fig_3_3():
    # Figure with expected maturities from 1 to 50
    from ER1 import prem as prem1
    from ER2 import prem as prem2
    from ER5 import prem as prem5
    data1 = np.array(prem1)
    data2 = np.array(prem2)
    data5 = np.array(prem5)
    E_mat = data1[:,0]
    prem1 = data1[:,1]
    prem2 = data2[:,1]
    prem5 = data5[:,1]

    fig, ax = plt.subplots()
    ax.plot(E_mat, prem1, label="r(0) = -1%")
    ax.plot(E_mat, prem2, label="r(0) = -2.5%")
    ax.plot(E_mat, prem5, label="r(0) = -5")
    plt.xlabel('Expected maturity (1 / $\lambda$)')
    plt.ylabel('Premium (%)')
    plt.legend()
    plt.show()

    # Figure with expected maturities from 1 to 500
    from ER12 import prem as prem12
    from ER22 import prem as prem22
    from ER52 import prem as prem52
    data12 = np.array(prem12)
    data22 = np.array(prem22)
    data52 = np.array(prem52)
    E_mat2 = data12[:,0]
    prem12 = data12[:,1]
    prem22 = data22[:,1]
    prem52 = data52[:,1]

    fig, ax = plt.subplots()
    ax.plot(E_mat2, prem12, label="r(0) = -1%")
    ax.plot(E_mat2, prem22, label="r(0) = -2.5%")
    ax.plot(E_mat2, prem52, label="r(0) = -5")
    plt.xlabel('Expected maturity (1 / $\lambda$)')
    plt.ylabel('Premium (%)')
    plt.legend()
    plt.show()


# Fig 4.1: Cash in Fixed supply versus long-dated Treasury STRIPs
def fig_4_1():
    # Prep data
    from fig41 import res
    data = np.array(res)
    nominalShortRate = data[0]
    ZCBPrice = data[1]
    cashPrice = data[2]
    cashPriceER = data[3]

    # Plot
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, ZCBPrice, label="30-year Treasury strip")
    ax.plot(nominalShortRate*100, cashPrice, label="Cash in fixed supply")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Price')
    plt.title("No early redemption")
    plt.legend()
    plt.show()

    # Plot early redemption
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, ZCBPrice, label="30-year Treasury strip")
    ax.plot(nominalShortRate*100, cashPriceER, label="Cash in fixed supply")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Price')
    plt.title("Early redemption ($\lambda$ = 2%)")
    plt.legend()
    plt.show()


# Fig 4.2: Duration of cash in xed supply.
def fig_4_2():
    # Prep data
    from fig42 import durs
    data = np.array(durs)
    nominalShortRate = data[0]
    ZCBDur = data[1]
    cashDur = data[2]
    cashDurER = data[3]
    print(ZCBDur)
    print(cashDur)
    print(cashDurER)
    # Plot
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, ZCBDur, label="30-year Treasury STRIP")
    ax.plot(nominalShortRate*100, cashDur, label="Cash in fixed supply")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Duration')
    plt.title("No early redemption")
    plt.legend()
    plt.show()

    # Plot early redemption
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, ZCBDur, label="30-year Treasury STRIP")
    ax.plot(nominalShortRate*100, cashDurER, label="Cash in fixed supply")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Duration')
    plt.title("Early redemption ($\lambda$ = 2%)")
    plt.legend()
    plt.show()

def fig_8():
    from out8 import bs
    data = np.array(bs)
    mat = data[:,0]
    bs = data[:,1]

    fig, ax = plt.subplots()
    ax.plot(mat, bs, label="r(0) = 0")
    plt.xlabel('T')
    plt.ylabel('barrier')
    plt.legend()
    plt.show()
