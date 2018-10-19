import random
import pandas


def get_data_random(filename):
    df = pandas.read_csv(filename)
    df2 = df.set_index("CLASS")
    data = df2.loc["toxic", "DNASEQ"]
    # grab some random samples
    sample = random.sample(range(len(data)), 400)
    dna = []
    for rand in sample:
        dna.append(data.iloc[rand])
    return dna


dna = get_data_random('peptidesWithDNA.csv')

