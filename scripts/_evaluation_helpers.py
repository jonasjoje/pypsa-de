import os
import pypsa

def load_networks_from_path_list(path_list):
    nn = {}
    for filepath in path_list:
        parts = filepath.split(os.sep)
        run = parts[parts.index("networks") - 1]
        year = int(filepath.split("_")[-1].split(".")[0])
        net = pypsa.Network(filepath)
        nn.setdefault(run, {})[year] = net
    return nn