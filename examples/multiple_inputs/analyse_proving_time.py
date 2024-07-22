import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


if __name__ == "__main__":
    mova_data = pd.read_csv('mova_point_prove_times.csv')
    powers = mova_data['pow'].tolist()
    prove_times_mova = mova_data['prove_time'].tolist()

    nova_data = pd.read_csv('nova_prove_times.csv')
    prove_times_nova = nova_data['prove_time'].tolist()





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
    
