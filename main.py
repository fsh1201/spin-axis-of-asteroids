from matplotlib import pyplot as plt

for k in range(0,320,20):
    plt.figure()
    for n in range(k, k+20):
        plt.subplot(4,5,n%20+1)
        b = str(n)
        a = 'E:\\dcf_test\\'
        c = '1.txt'
        a1 = a + b + c
        f = open(a1)
        data = f.read().splitlines()
        f = len(data)
        time = []
        bright = []
        for i in range(f):
            m = data[i].split()
            time.append(float(m[0]))
            bright.append(float(m[1]))

        plt.scatter(time, bright, c='r')

        c = '2.txt'
        a2 = a + b + c

        f = open(a2)
        data = f.read().splitlines()
        f = len(data)
        time = []
        bright = []
        for i in range(f):
            m = data[i].split()
            time.append(float(m[0]))
            bright.append(float(m[1]))

        plt.plot(time, bright, c='b')

    plt.show()
