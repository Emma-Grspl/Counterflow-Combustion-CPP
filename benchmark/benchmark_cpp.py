from pathlib import Path
import csv
import statistics
import subprocess
import time


EXECUTABLE = Path(
    "build/counterflow_combustion"
)

STEPS = 2999
WARMUP_RUNS = 1
MEASURED_RUNS = 7

OUTPUT = Path(
    "benchmark/cpp_timings.csv"
)


def run_once():
    start = time.perf_counter()

    subprocess.run(
        [
            str(EXECUTABLE),
            str(STEPS),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    end = time.perf_counter()

    return end - start


def main():
    if not EXECUTABLE.exists():
        raise SystemExit(
            f"Missing executable: {EXECUTABLE}"
        )

    print("C++ solver benchmark")
    print("====================")
    print(
        f"Warm-up runs  : {WARMUP_RUNS}"
    )
    print(
        f"Measured runs : {MEASURED_RUNS}"
    )
    print()

    for i in range(WARMUP_RUNS):
        elapsed = run_once()

        print(
            f"warm-up {i + 1}: "
            f"{elapsed:.6f} s"
        )

    timings = []

    print()

    for i in range(MEASURED_RUNS):
        elapsed = run_once()

        timings.append(elapsed)

        print(
            f"run {i + 1:2d}: "
            f"{elapsed:.6f} s"
        )

    mean = statistics.mean(
        timings
    )

    median = statistics.median(
        timings
    )

    std = (
        statistics.stdev(timings)
        if len(timings) > 1
        else 0.0
    )

    minimum = min(timings)
    maximum = max(timings)

    print()
    print("Summary")
    print("-------")

    print(
        f"mean   : {mean:.6f} s"
    )

    print(
        f"median : {median:.6f} s"
    )

    print(
        f"std    : {std:.6f} s"
    )

    print(
        f"min    : {minimum:.6f} s"
    )

    print(
        f"max    : {maximum:.6f} s"
    )

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with OUTPUT.open(
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)

        writer.writerow([
            "run",
            "seconds",
        ])

        for i, elapsed in enumerate(
            timings,
            start=1,
        ):
            writer.writerow([
                i,
                elapsed,
            ])

    print()
    print(
        f"Raw timings: {OUTPUT}"
    )


if __name__ == "__main__":
    main()
