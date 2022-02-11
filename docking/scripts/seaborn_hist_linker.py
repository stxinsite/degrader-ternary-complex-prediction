import os
import sys

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def main(filename):

    data = pd.read_csv(filename, sep='\t', lineterminator='\n')

    sns.set(font_scale=1.3)

    h1 = sns.histplot(x=data.length, kde=True, color = "red",)
    h1.text(4.5, 120, "Linker 651: linker\nmean = 5.0 $\AA$ \nsd = 0.8 $\AA$", fontsize=25)
    h1.set_xlabel("end-to-end distance [$\AA$]", fontsize=25)
    h1.set_ylabel("Count", fontsize=25)

    plt.tight_layout()
    plt.show()


if "__main__" == __name__:
    filename = sys.argv[1]
    main(filename)
