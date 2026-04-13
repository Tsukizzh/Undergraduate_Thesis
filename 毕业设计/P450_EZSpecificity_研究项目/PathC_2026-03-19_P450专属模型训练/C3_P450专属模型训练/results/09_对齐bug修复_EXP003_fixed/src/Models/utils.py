import sys

root_dir = "/work/yufeng/2022/enzyme_specificity"

sys.path.append(f"{root_dir}/src")

def load_model(config):
    from Models.dlkcat import DLKcat
    from Models.cpi import CPI
    from Models.ss import SS

    if config.model.name == 'DLKcat':
        return DLKcat(config)
    elif config.model.name == 'cpi':
        return CPI(config)
    elif config.model.name == 'structure_sequence':
        return SS(config)
    else:
        raise ValueError(f'Unknown model name: {config.model.name}')

def increase_count(x, a, cnt = 1):
    if a in x:
        x[a] += cnt
    else:
        x[a] = cnt

def get_count(x, a):
    if a in x:
        return x[a]
    else:
        return 0

def create_fan_chart(datas, save_dir):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    plt.clf()
    data = datas
    labels = np.arange(len(datas))
    colors = sns.color_palette('pastel')
    plt.pie(data, labels=labels,colors = colors, autopct = '%0.01f%%')
    plt.savefig(save_dir)

def exists(val):
    return val is not None

def default(val, d):
    return val if exists(val) else d

def stable_softmax(t, dim = -1):
    t = t - t.amax(dim = dim, keepdim = True)
    return t.softmax(dim = dim)