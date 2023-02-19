import time
from typing import List

import tqdm
from assess_prealloc_list import MyItem, generate_random_items


def main_b(num_items: int) -> List[MyItem]:
    retl: List[MyItem] = [MyItem.blank() for _ in range(num_items)]
    for i, item in enumerate(generate_random_items(num_items)):
        retl[i] = item
    return retl


def main_bn(num_items: int) -> List[MyItem]:
    retl: List[MyItem] = [MyItem.blank_none() for _ in range(num_items)]
    for i, item in enumerate(generate_random_items(num_items)):
        retl[i] = item
    return retl


def main(num_items: int) -> List[MyItem]:
    retl: List[MyItem] = list(generate_random_items(num_items))
    return retl


if __name__ == "__main__":
    with open("out.csv", "w") as writer:
        writer.write(",".join(("FN", "POW", "REP", "TIME")) + "\n")
        for pow in range(1, 7):
            for replicate in tqdm.tqdm(range(100), desc=str(pow)):
                for fn in main, main_b, main_bn:
                    ts = time.time()
                    fn(10 ** pow)
                    te = time.time()
                    writer.write(",".join((
                        fn.__name__,
                        str(pow),
                        str(replicate),
                        str(te - ts)
                    )) + "\n")
