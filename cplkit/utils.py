"""Small utility functions for command-line reporting."""

from __future__ import annotations

import math
import sys
import time


class Timer:
    def __init__(self) -> None:
        self.t0 = time.perf_counter()

    def elapsed(self) -> float:
        return time.perf_counter() - self.t0

    @staticmethod
    def fmt(seconds: float) -> str:
        if seconds < 60:
            return f"{seconds:5.1f}s"
        if seconds < 3600:
            m = int(seconds // 60)
            s = seconds - 60 * m
            return f"{m:d}m{s:04.1f}s"
        h = int(seconds // 3600)
        m = int((seconds - 3600 * h) // 60)
        s = seconds - 3600 * h - 60 * m
        return f"{h:d}h{m:02d}m{s:04.1f}s"


def progress_line(prefix: str, i: int, n: int, t0: float) -> str:
    now = time.perf_counter()
    done = i / max(n, 1)
    elapsed = now - t0
    eta = (elapsed / done - elapsed) if done > 0 else float("inf")
    bar_w = 26
    fill = int(bar_w * done)
    bar = "[" + "#" * fill + "-" * (bar_w - fill) + "]"
    eta_s = "  ETA " + Timer.fmt(eta) if math.isfinite(eta) else ""
    return f"{prefix} {bar} {i}/{n}  elapsed {Timer.fmt(elapsed)}{eta_s}"


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)
