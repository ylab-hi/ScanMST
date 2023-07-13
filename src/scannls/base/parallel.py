"""Module for parallel worker.

@Filename:    parallel.py
@Author:      YangyangLi
@Time:        12/15/21 1:58 PM
"""
import multiprocessing
import os
from collections.abc import Callable
from concurrent import futures
from typing import Any

from scannls.type import LoggerType


class ParallelWorker:
    """ParallelWorker class is used to run function in parallel.

    args include the unique parameter of the function and  keyword arguments include
    the common parameters of the function


    :param func: the function to be run in parallel
    :param n_jobs: the number of jobs to run in parallel
    :param logger: the logger object

    :Example:

    >>> from loguru import logger
    >>> def func(x, *, y=1):
    ...     z = x + y
    ...     return z
    >>> args, kwargs = [1, 2, 3], {'y': 4}
    >>> n_jobs = 3
    >>> parallel_worker = ParallelWorker(func=func, logger=logger, n_jobs=n_jobs)
    >>> result = parallel_worker.run(*args, **kwargs)
    >>> result
    {1: 5, 2: 6, 3: 7}
    """

    def __init__(
        self,
        func: Callable[..., Any],
        logger: LoggerType,
        n_jobs: int = 1,
    ) -> None:
        """Initialize the ParallelWorker class."""
        self.func = func
        self.logger = logger

        self.n_jobs = self.setter_n_jobs(n_jobs)

    def setter_n_jobs(self, n_jobs: int) -> int:
        """Set the number of jobs to run in parallel in terms of cpu cores."""
        current_max_processor = os.cpu_count()

        if current_max_processor is None:
            current_max_processor = n_jobs

        if n_jobs > current_max_processor:
            self.logger.warning(
                f"ParallelWorker: {n_jobs} > current_max_processor {current_max_processor}",
            )
            return current_max_processor

        return n_jobs  # the max processor is decided by ProcessPoolExecutor

    def run(self, *args: Any, **kwargs: Any) -> dict[str, Any]:
        """Using concurrent.future to parallel process."""
        tasks = {}
        result = {}
        self.logger.info(f"ParallelWorker: {self.n_jobs} jobs")
        m = multiprocessing.Manager()
        lock = m.Lock()  # add lock to protect blat log
        with futures.ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            for key in args:
                self.logger.debug(f"ParallelWorker: {key} submitted")
                future = executor.submit(self.func, key, lock, **kwargs)
                tasks[future] = key

            for future in futures.as_completed(tasks):
                self.logger.trace(f"ParallelWorker: {tasks[future]} done")
                key = tasks[future]
                result[key] = future.result()
        return result

    def map(self, *iterables, timeout=None, chunksize=1) -> Any:
        """Using concurrent.futures to parallel process."""
        with futures.ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            return executor.map(
                self.func,
                *iterables,
                timeout=timeout,
                chunksize=chunksize,
            )
