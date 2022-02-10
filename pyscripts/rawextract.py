import os
dataset = "artificial_human_chr22_DeepSimu"

path = "/Users/aaronphilip/ScienceFair2022/DeepSimulator/DeepSimulator/%s/fast5" % dataset

ext = ('.fast5')
files = os.scandir(path)

for file in files:
    if file.path.endswith(ext):
        os.system("h5ls -d ./%s/fast5/%s/Raw/Reads/Read_981/Signal >> ./signaldata/%s_signals.py" % (dataset,file.name, dataset))

