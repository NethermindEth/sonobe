import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Later need to read this from csv
    powers = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    prove_times_mova = [1.091167, 1.415209, 1.959292, 3.783584, 6.825875, 12.65875, 27.475917, 49.278416, 104.126292, 210.376334, 423.0625, 802.578584, 1817.999292, 5518.4255]
    prove_times_nova = [3.464959, 4.878416, 7.162625, 12.449458, 25.237625, 45.0025, 71.196208, 170.822917, 284.674875, 622.635709, 1147.108875, 2197.209292, 4823.36525, 11638.812291]

    powers = list(map(str, powers))
    powers = ["2^" + pw for pw in powers]

    ratios = np.array(prove_times_nova) / np.array(prove_times_mova)
    print(ratios.mean())

    # Proving time in log scale
    fig = plt.figure(figsize =(10, 7))
    plt.yscale("log")
    x_axis = np.arange(len(powers)) 
    plt.bar(x_axis-0.2, prove_times_mova, 0.4, label="Mova")
    plt.bar(x_axis+0.2, prove_times_nova, 0.4, label="Nova")

    plt.xticks(x_axis, powers) 
    plt.title("Analyzing proving times")
    plt.xlabel("Input size")
    plt.ylabel("Proving time (ms) in log scale")

    plt.savefig("prove_time_log.png")
    plt.legend()
    plt.show()

    # Proving time absolute
    fig = plt.figure(figsize =(10, 7))

    plt.bar(x_axis-0.2, prove_times_mova, 0.4, label="Mova")
    plt.bar(x_axis+0.2, prove_times_nova, 0.4, label="Nova")

    plt.xticks(x_axis, powers) 
    plt.title("Analyzing proving times")
    plt.xlabel("Input size")
    plt.ylabel("Proving time (ms)")
    
    plt.savefig("prove_time.png")
    plt.legend()
    plt.show()
    
