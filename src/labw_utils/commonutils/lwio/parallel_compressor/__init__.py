"""
NOT FINISHED -- DO NOT USE
"""

import queue
import threading

from labw_utils.typing_importer import Callable, BinaryIO, Optional

COMPRESS_FUNC_TYPE = Callable[[bytes], bytes]


class ParallelCompressorThreadType(threading.Thread):
    _in_data: bytes
    _out_data: Optional[bytes]
    _compress_func: COMPRESS_FUNC_TYPE

    def __init__(self, in_data: bytes, compress_func: COMPRESS_FUNC_TYPE):
        super().__init__()
        self._in_data = in_data
        self._out_data = None
        self._compress_func = compress_func

    def run(self):
        self._out_data = self._compress_func(self._in_data)

    @property
    def out_data(self) -> bytes:
        if self._out_data is None:
            raise RuntimeError("Data not ready!")
        return self._out_data


class ParallelReader(threading.Thread):
    _src_io: BinaryIO
    _out_queue: queue.Queue
    _chunksize: int
    _compress_func: COMPRESS_FUNC_TYPE

    def __init__(self, src_io: BinaryIO, out_queue: queue.Queue, chunksize: int, compress_func: COMPRESS_FUNC_TYPE):
        super().__init__()
        self._src_io = src_io
        self._out_queue = out_queue
        self._chunksize = chunksize
        self._compress_func = compress_func

    def run(self):
        print("Reader START")
        while True:
            chunk = self._src_io.read(self._chunksize)
            if not chunk:
                break
            compressor = ParallelCompressorThreadType(in_data=chunk, compress_func=self._compress_func)
            compressor.start()
            self._out_queue.put(compressor)
        print("Reader FIN")


class _ParallelCompressorWriter(threading.Thread):
    _dst_io: BinaryIO
    _in_queue: queue.Queue
    _is_finished: bool

    def __init__(self, dst_io: BinaryIO, in_queue: queue.Queue):
        super().__init__()
        self._dst_io = dst_io
        self._in_queue = in_queue
        self._is_finished = False

    def run(self):
        print("Writer START")
        while True:
            try:
                compressor: ParallelCompressorThreadType = self._in_queue.get_nowait()
            except queue.Empty:
                if self._is_finished:
                    break
                continue
            compressor.join()
            self._dst_io.write(compressor.out_data)
        print("Writer FIN")

    def terminate(self):
        self._is_finished = True


def parallel_compress(
    src_io: BinaryIO, dst_io: BinaryIO, compress_func: COMPRESS_FUNC_TYPE, chunksize: int = 4 * 1024 * 1024
):
    q = queue.Queue()
    reader_thread = ParallelReader(src_io, q, chunksize, compress_func)
    reader_thread.start()
    writer_thread = _ParallelCompressorWriter(dst_io, q)
    writer_thread.start()
    reader_thread.join()
    writer_thread.terminate()
    writer_thread.join()
