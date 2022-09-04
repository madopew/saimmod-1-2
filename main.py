import math
import random

import matplotlib.pyplot as plt


def lehmerStep(r, a, m):
    return (a * r) % m


def lehmer(r, a, m, n):
    if a >= m:
        raise ValueError("a must be less than m")

    result = []
    for i in range(n):
        r = lehmerStep(r, a, m)
        result.append(r / m)
    return result


def rndLehmer(n, nue):
    ls = []
    for i in range(nue):
        r = random.randint(2, 100)
        a = random.randint(1000, 5000)
        m = random.randint(10000000, 20000000)
        ls.append(lehmer(r, a, m, n))
    return ls


def uniform(a, b, n):
    return [a + (b - a) * x for x in rndLehmer(n, 1)[0]]


def normal(avg, sigma, n, nue=6):
    ls = rndLehmer(n, nue)

    x = []
    for i in range(n):
        xElement = 0
        for j in range(nue):
            xElement += ls[j][i]
        sqrt = (12 / nue) ** 0.5
        x.append(avg + sigma * sqrt * (xElement - nue / 2))

    return x


def exp(lambd, n):
    return [-math.log(1 - x) / lambd for x in rndLehmer(n, 1)[0]]


def gamma(lambd, nue, n):
    ls = rndLehmer(n, nue)

    x = []
    for i in range(n):
        xElement = 1
        for j in range(nue):
            xElement *= ls[j][i]
        x.append(-math.log(xElement) / lambd)

    return x


def triangle(a, b, n, tp="first"):
    ls = rndLehmer(n, 2)

    x = []
    for i in range(n):
        if tp == "first":
            x.append(a + (b - a) * max(ls[0][i], ls[1][i]))
        elif tp == "second":
            x.append(a + (b - a) * min(ls[0][i], ls[1][i]))
        else:
            raise ValueError("unknown type")

    return x


def simpson(a, b, n):
    uni1 = uniform(a / 2, b / 2, n)
    uni2 = uniform(a / 2, b / 2, n)
    return [x + y for x, y in zip(uni1, uni2)]


def mathProperties(values):
    average = sum(values) / len(values)
    variance = sum([(x - average) ** 2 for x in values]) / len(values)
    sigma = variance ** 0.5
    return average, variance, sigma


def sequenceProperties(values, eps=1e-6):
    anchor = values[-1]
    periodStart = -1
    periodEnd = -1

    for i in range(len(values) - 1):
        if abs(anchor - values[i]) < eps:
            if periodStart == -1:
                periodStart = i
            else:
                periodEnd = i
                break

    if periodStart == -1:
        raise ValueError("period start not found")
    if periodEnd == -1:
        raise ValueError("period end not found")

    period = periodEnd - periodStart

    aperiodicLength = 0
    for i in range(periodStart):
        if abs(values[i] - values[i + period]) < eps:
            break
        aperiodicLength += 1

    k = 0
    for i in range(0, len(values), 2):
        if (values[i] ** 2 + values[i + 1] ** 2) <= 1:
            k += 1

    return aperiodicLength, period, k


def displayProperties(values, name):
    average, variance, sigma = mathProperties(values)
    print("==================== {} ====================".format(name))
    print("Average: {}".format(average))
    print("Variance: {}".format(variance))
    print("Sigma: {}".format(sigma))
    print("============================================")


def lab1():
    # good ones
    # r = 15
    # a = 1643
    # m = 12031278
    # n = 1_000_000

    r = 15
    a = 1643
    m = 12031278
    n = 1_000_000

    values = lehmer(r, a, m, n)
    aperiodicLength, period, k = sequenceProperties(values)

    print("1/12 = {}".format(1 / 12))
    print("1/sqrt(12) = {}".format(1 / 12 ** 0.5))
    print("pi/4 = {}".format(3.141592653589793 / 4))
    print("")
    displayProperties(values, "Lehmer")
    print("Aperiodic length: {}".format(aperiodicLength))
    print("Period: {}".format(period))
    print("2K/N: {}".format(2 * k / n))

    plt.hist(values, bins=20)
    plt.show()


def lab2():
    n = 100_000

    uniformSeq = uniform(2, 5, n)
    displayProperties(uniformSeq, "Uniform")
    normalSeq = normal(3, 1, n)
    displayProperties(normalSeq, "Normal")
    expSeq = exp(2, n)
    displayProperties(expSeq, "Exp")
    gammaSeq = gamma(1, 6, n)
    displayProperties(gammaSeq, "Gamma")
    triangleSeq = triangle(-4, 5, n, "second")
    displayProperties(triangleSeq, "Triangle")
    simpsonSeq = simpson(4, 8, n)
    displayProperties(simpsonSeq, "Simpson")

    plt.hist(uniformSeq, bins=20)
    plt.show()
    plt.hist(normalSeq, bins=20)
    plt.show()
    plt.hist(expSeq, bins=20)
    plt.show()
    plt.hist(gammaSeq, bins=20)
    plt.show()
    plt.hist(triangleSeq, bins=20)
    plt.show()
    plt.hist(simpsonSeq, bins=20)
    plt.show()


def main():
    lab1()
    lab2()


if __name__ == "__main__":
    main()
