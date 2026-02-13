import os
dirname = ("\\").join(os.path.dirname(__file__).split("\\")[:-1])
filename = os.path.join(dirname, 'files\\test.tsv')
print(filename)
