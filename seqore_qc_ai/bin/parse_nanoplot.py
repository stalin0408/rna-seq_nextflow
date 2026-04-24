#!/usr/bin/env python3

import re
import json
import sys

file = sys.argv[1]

def clean(x):
    return float(x.replace(",", "").strip())

metrics = {}

with open(file) as f:
    for line in f:
        line = line.strip()

        # Core metrics
        if line.startswith("Mean read length"):
            metrics["mean_read_length"] = clean(line.split(":")[1])

        elif line.startswith("Mean read quality"):
            metrics["mean_quality"] = clean(line.split(":")[1])

        elif line.startswith("Median read length"):
            metrics["median_read_length"] = clean(line.split(":")[1])

        elif line.startswith("Median read quality"):
            metrics["median_quality"] = clean(line.split(":")[1])

        elif line.startswith("Number of reads"):
            metrics["num_reads"] = int(clean(line.split(":")[1]))

        elif line.startswith("Read length N50"):
            metrics["n50"] = clean(line.split(":")[1])

        elif line.startswith("STDEV read length"):
            metrics["stdev_length"] = clean(line.split(":")[1])

        elif line.startswith("Total bases"):
            metrics["total_bases"] = clean(line.split(":")[1])

        # Quality thresholds (>Q10, >Q15, etc.)
        elif line.startswith(">Q"):
            match = re.match(r">Q(\d+):\s+([\d,]+)\s+\(([\d\.]+)%\)", line)
            if match:
                q = match.group(1)
                count = int(match.group(2).replace(",", ""))
                pct = float(match.group(3))

                metrics[f"reads_above_Q{q}"] = {
                    "count": count,
                    "percent": pct
                }

# Output JSON
print(json.dumps(metrics, indent=2))