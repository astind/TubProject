import pandas
import random


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


def print_output(best_ans):
    out_string = ""
    for ans in best_ans:
        print(ans)
        out_string += ans + '\n'
    out_string = out_string.strip()
    out_file = open("randomData.txt", "w")
    out_file.write(out_string)
    out_file.close()


data = get_data_random("peptidesWithDNA.csv")
print_output(data)


