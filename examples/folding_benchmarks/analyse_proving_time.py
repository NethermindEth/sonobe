import pandas as pd


if __name__ == "__main__":
    mova_data = pd.read_csv('mova_prove_times.csv')
    nova_data = pd.read_csv('nova_prove_times.csv')
    hp_data = pd.read_csv('hypernova_prove_times.csv')

    mova_prove_time = mova_data.groupby('pow')['prove_time'].mean() / 1000
    nova_prove_time = nova_data.groupby('pow')['prove_time'].mean() / 1000
    hp_prove_time = hp_data.groupby('pow')['prove_time'].mean() / 1000

    print("Mova", mova_prove_time)
    print("Nova", nova_prove_time)
    print("Hypernova", hp_prove_time)






    
