from scipy.stats import chisquare

observed = [100e-20, 200e-20, 300e-20, 400e-20, 500e-20]
expected = [3e-20, 3e-20, 3e-20, 3e-20, 3e-20]

chisquare = chisquare(observed, expected)

print(chisquare, pvalue)