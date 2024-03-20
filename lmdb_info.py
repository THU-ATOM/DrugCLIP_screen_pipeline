import lmdb
import pickle
from tqdm import tqdm
from random import shuffle
def read_lmdb(lmdb_path):
    env = lmdb.open(
        lmdb_path,
        subdir=False,
        readonly=True,
        lock=False,
        readahead=False,
        meminit=False,
        max_readers=256,
    )
    txn = env.begin()
    keys = list(txn.cursor().iternext(values=False))
    set_count = {}
    for idx in tqdm(keys):
        datapoint= txn.get(idx)
        datapoint_pickled = pickle.loads(datapoint)
        try:
            set_count[datapoint_pickled['subset']] += 1
        except:
            set_count.update({datapoint_pickled['subset']:1})
    env.close()
    #sort set_count by key name
    set_count = dict(sorted(set_count.items()))
    return set_count

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('lmdb_path', type=str, help='lmdb path')
    args = parser.parse_args()
    set_count = read_lmdb(args.lmdb_path)
    print(set_count)