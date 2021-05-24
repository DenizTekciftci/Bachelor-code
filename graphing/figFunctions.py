import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import seaborn as sns

w = 7
h = 5

# Fig P
def fig_P():
    from data.out21 import zcbs as P1
    from data.out22 import zcbs as P2
    from data.out212 import times as t1
    from data.out222 import times as t2
    # Regular sampling data
    data1 = np.array(P1)
    mean1 = data1[:,0]
    std1 = data1[:,1]
    loConf1 = data1[:,2]
    hiConf1 = data1[:,3]

    # Antithetic sampling data
    data2 = np.array(P2)
    mean2 = data2[:,0]
    std2 = data2[:,1]
    loConf2 = data2[:,2]
    hiConf2 = data2[:,3]

    maturities = range(51)

    # Plot regular sampling P(0,T)
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(maturities, mean1, color="orange")
    ax.fill_between(maturities, loConf1, hiConf1, alpha=0.8)
    plt.xlabel('Maturities')
    plt.ylabel('Price of Zero Coupon Bond')
    plt.grid(axis="both")
    plt.show()

    # Plot AT sampling P(0,T)
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(maturities, mean2, color="orange")
    ax.fill_between(maturities, loConf2, hiConf2, alpha=0.8)
    plt.xlabel('Maturities')
    plt.ylabel('Price of Zero Coupon Bond')
    plt.grid(axis="both")
    plt.show()

    # Put into np arrays so we can divide easily
    t1 = np.array(t1)
    t2 = np.array(t2)
    #Plot efficiency ratios and time ratios
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(maturities, (std2/std1)**2, color="orange", label="Efficiency ratio")
    ax.plot(maturities, (t1/t2), color="blue", label="Time ratio")
    plt.xlabel('Maturities')
    plt.ylabel('Ratio')
    plt.grid(axis="both")
    plt.legend()
    plt.show()


def fig_3_1():
    from data.fig31 import volSurface as data
    data = np.array(data)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    for i in range(50, 70):
        print(len(x)/i)

    x = np.reshape(x, (60, 30))
    y = np.reshape(y, (60, 30))
    z = np.reshape(z, (60, 30))



    print(y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('Nominal short rate (%)')
    ax.set_ylabel('Short rate volatility (%)')
    ax.set_zlabel('Premium (%)')
    ax.plot_surface(x*100,y*100,z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    plt.show()

# Fig 3.3: Term structure of cash
def fig_3_3():
    # Figure with expected maturities from 1 to 50
    from data.ER1 import prem as prem1
    from data.ER2 import prem as prem2
    from data.ER5 import prem as prem5
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
    plt.grid(axis="both")
    plt.show()

    # Figure with expected maturities from 1 to 500
    from data.ER12 import prem as prem12
    from data.ER22 import prem as prem22
    from data.ER52 import prem as prem52
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
    plt.grid(axis="both")
    plt.show()




# Fig 4.1: Cash in Fixed supply versus long-dated Treasury STRIPs
def fig_4_1():
    # Prep data
    from data.fig41 import res
    data = np.array(res)
    nominalShortRate = data[0]
    ZCBPrice = data[1]
    cashPrice = data[2]
    cashPriceER1 = data[3]
    cashPriceER2 = data[4]
    cashPriceER3 = data[5]
    cashPriceER4 = data[6]

    # Plot early redemption
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, ZCBPrice, label="30-year Treasury strip")
    ax.plot(nominalShortRate*100, cashPrice, label="Cash ($\lambda$ = 0%)")
    ax.plot(nominalShortRate*100, cashPriceER1, label="Cash ($\lambda$ = 2%)")
    ax.plot(nominalShortRate*100, cashPriceER2, label="Cash ($\lambda$ = 4%)")
    ax.plot(nominalShortRate*100, cashPriceER3, label="Cash ($\lambda$ = 8%)")
    ax.plot(nominalShortRate*100, cashPriceER4, label="Cash ($\lambda$ = 16%)")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Price')
    plt.title("Cash in fixed supply")
    plt.legend()
    plt.grid(axis="both")
    plt.show()


# Fig 4.2: Duration of cash in xed supply.
def fig_4_2():
    # Prep data
    from data.fig42 import durs
    data = np.array(durs)
    nominalShortRate = data[0]
    stripDur = data[1]
    cashDur = data[2]
    cashDurER1 = data[3]
    cashDurER2 = data[4]
    cashDurER3 = data[5]
    cashDurER4 = data[6]

    # # Plot
    # fig, ax = plt.subplots()
    # ax.plot(nominalShortRate*100, stripDur, label="30-year Treasury STRIP")
    # plt.xlabel('Nominal short rate (%)')
    # plt.ylabel('Duration')
    # plt.title("No early redemption")
    # plt.legend()
    # plt.show()

    # Plot early redemption
    fig, ax = plt.subplots()
    ax.plot(nominalShortRate*100, stripDur, label="30-year Treasury STRIP")
    ax.plot(nominalShortRate*100, cashDur, label="Cash ($\lambda$ = 0%)")
    ax.plot(nominalShortRate*100, cashDurER1, label="Cash ($\lambda$ = 2%)")
    # ax.plot(nominalShortRate*100, cashDurER2, label="Cash ($\lambda$ = 4%)")
    # ax.plot(nominalShortRate*100, cashDurER3, label="Cash ($\lambda$ = 8%)")
    # ax.plot(nominalShortRate*100, cashDurER4, label="Cash ($\lambda$ = 16%)")
    plt.xlabel('Nominal short rate (%)')
    plt.ylabel('Duration')
    plt.title("Cash in fixed supply")
    plt.legend()
    plt.grid(axis="both")
    plt.show()

# Draw paths of r
def fig_5_2():
    from data.out52 import sims as data
    data = np.array(data)

    w = 9
    h = 7
    fig, ax = plt.subplots(figsize=(w, h))
    plt.xlabel('Time')
    plt.ylabel('Nominal short rate (%)')
    for i in range(len(data)//2):
        ax.plot(data[i*2+1], data[i*2]*100)

    plt.grid(axis="both")
    plt.show()

# Draw paths of r
def fig_5_3():
    from data.out53 import sims as data
    data = np.array(data)

    w = 9
    h = 7
    fig, ax = plt.subplots(figsize=(w, h))
    plt.xlabel('Time')
    plt.ylabel('Nominal short rate (%)')
    for i in range(len(data)//2):
        ax.plot(data[i*2+1], data[i*2]*100)
        
    plt.grid(axis="both")
    plt.show()


def fig_8():
    from data.out8 import B
    data = np.array(B)

    # Removes very last element of T = 5 and T = 50, which are wrong
    del data[0][-1]
    del data[1][-1]
    del data[8][-1]
    del data[9][-1]

    mat1 = np.array(data[0])
    bs1 = np.array(data[1])
    mat2 = np.array(data[2])
    bs2 = np.array(data[3])
    mat3 = np.array(data[4])
    bs3 = np.array(data[5])
    mat4 = np.array(data[6])
    bs4 = np.array(data[7])
    mat5 = np.array(data[8])
    bs5 = np.array(data[9])


    fig, ax = plt.subplots()
    ax.plot(mat1, bs1 * 100, label="T = 5")
    ax.plot(mat2, bs2 * 100, label="T = 10")
    ax.plot(mat3, bs3 * 100, label="T = 15")
    ax.plot(mat4, bs4 * 100, label="T = 30")
    ax.plot(mat5, bs5 * 100, label="T = 50")
    plt.xlabel('Time')
    plt.ylabel('Barrier (%)')
    plt.legend()
    plt.grid(axis="both")
    plt.show()

# Plot time to barrier
def fig_6():
    from data.out6 import taus
    data = np.array(taus)
    mean = np.mean(data)
    sns.distplot(data)
    print(sns.distplot(data))
    plt.axvline(x = mean, linestyle='--', color='blue', label = "mean = {}".format(round(mean,2)))
    plt.xlabel('Time')
    plt.ylabel('Relative frequncy')
    plt.legend()
    plt.grid(axis="both")
    plt.show()


def fig_8_1():
    from data.out81 import C
    fig, ax = plt.subplots()
    ax.plot(range(5, 105, 5), np.array(C), label="$C^T(0)$")
    ax.plot(range(5, 105, 5), [1.44195]*20, label="$C^\infty(0)$")
    plt.xlabel('Time')
    plt.ylabel('ZCB-Price')
    plt.grid(axis="both")
    plt.legend()
    plt.show()



def fig_8_2():
    from data.out82 import B
    data = np.array(B)
    data = np.array(data)
    sigma = np.array(data[:,1])
    b = np.array(data[:,0])


    print(b)

    fig, ax = plt.subplots()
    ax.plot(sigma * 100, b * 100)
    plt.xlabel('Short rate volatility (%)')
    plt.ylabel('Barrier (%)')
    plt.grid(axis="both")
    plt.show()


def fig_8_3():
    from data.out83 import C
    fig, ax = plt.subplots()
    ax.plot(range(10, 300, 20), np.array(C) , label="$C^T(0)$")
    ax.plot(range(10, 300, 20), [1.44195]*15, label="$C^\infty(0)$")
    plt.xlabel('Time')
    plt.ylabel('ZCB-Price')
    plt.grid(axis="both")
    plt.legend()
    plt.show()