import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("res_run-20250825T145232.csv")  # change the name
plt.plot(df["Vc"], df["energia"], "-o", markersize=2)
plt.xlabel("Vc")
plt.ylabel("calculated Energy")
plt.plot([min(df["Vc"]), max(df["Vc"])], [-1.990390] * 2, "--")
# Change constant line depending if J=3, J=1
# ground -2.467000 ; excited -1.990390
plt.show()
