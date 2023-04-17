import os.path
import time
from labw_utils.typing_importer import Type

import numpy as np
import tqdm

from naive_interval_engine import IntervalEngineType
from naive_interval_engine.ne_impl import NumExprIntervalEngine
from naive_interval_engine.np_impl import NumpyIntervalEngine
from naive_interval_engine.pd_impl import PandasIntervalEngine

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_FILE_DIR = os.path.join(
    FILE_DIR,
    "test",
    "test_naive_interval_engine"
)


def assess_nie(
        nie_type: Type[IntervalEngineType],
        gtf_regions_filepath: str,
        bam_regions_filepath: str,
        show_tqdm: bool
):
    with open(os.path.join(FILE_DIR, "log.log"), "a") as log_appender:
        def write_log(contents):
            log_appender.write(
                f"{time.asctime()}: "
                f"{nie_type.__name__} 4 "
                f"{os.path.basename(gtf_regions_filepath)} + "
                f"{os.path.basename(bam_regions_filepath)}: "
                f"{contents}\n"
            )
            log_appender.flush()

        write_log("start")
        timestamps = [time.time()]
        gtf_nie = nie_type(gtf_regions_filepath, show_tqdm=show_tqdm)
        bam_nie = nie_type(bam_regions_filepath, show_tqdm=show_tqdm)
        write_log("parsed")
        timestamps.append(time.time())
        _ = list(gtf_nie.matches(bam_nie, show_tqdm=show_tqdm))
        write_log("match1")
        timestamps.append(time.time())
        _ = list(bam_nie.matches(gtf_nie, show_tqdm=show_tqdm))
        write_log("match2")
        timestamps.append(time.time())
        _ = list(gtf_nie.overlaps(bam_nie, show_tqdm=show_tqdm))
        write_log("overlap1")
        timestamps.append(time.time())
        _ = list(bam_nie.overlaps(gtf_nie, show_tqdm=show_tqdm))
        write_log("overlap2")
        timestamps.append(time.time())
        timestamps_arr = np.array(timestamps)
        reta = timestamps_arr[1:] - timestamps_arr[:-1]
        write_log("fin")
        return reta


def assess_using_test_data(
        nie_type: Type[IntervalEngineType],
        n_trials: int = 100
):
    final_assess = np.zeros((n_trials, 5), dtype=float)
    for i in tqdm.tqdm(range(n_trials)):
        final_assess[i] = assess_nie(
            nie_type,
            os.path.join(TEST_FILE_DIR, "test_match.tsv"),
            os.path.join(TEST_FILE_DIR, "test_data.tsv"),
            show_tqdm=False
        )
    return np.sum(final_assess, axis=0)


def assess_using_real_data(nie_type: Type[IntervalEngineType]):
    return assess_nie(
        nie_type,
        os.path.join(FILE_DIR, "gtf_regions_chr1.tsv"),
        os.path.join(FILE_DIR, "bam_regions_chr1.tsv"),
        show_tqdm=True
    )


if __name__ == "__main__":
    with open(os.path.join(FILE_DIR, "result.tsv"), "w") as writer:
        writer.write("\t".join((
            "Implementation",
            "Data",
            "Create",
            "Match1",
            "Match2",
            "Overlap1",
            "Overlap2"
        )) + "\n")
        for nie_type in (
                NumpyIntervalEngine,
                PandasIntervalEngine,
                NumExprIntervalEngine,
                # IntervalTreeIntervalEngine, # Ultra slow
        ):
            test_data_result = assess_using_test_data(nie_type=nie_type)
            writer.write("\t".join((
                nie_type.__name__,
                "test",
                *map(str, test_data_result)
            )) + "\n")
            bulk_data_result = assess_using_real_data(nie_type=nie_type)
            writer.write("\t".join((
                nie_type.__name__,
                "real",
                *map(str, bulk_data_result)
            )) + "\n")
