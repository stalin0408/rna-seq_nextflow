#!/usr/bin/env python3

import json
import sys

before_file = sys.argv[1]
after_file = sys.argv[2]

with open(before_file) as f:
    before = json.load(f)

with open(after_file) as f:
    after = json.load(f)

def diff(a, b):
    if a is None or b is None:
        return None
    return round(b - a, 2)

report = {}

keys = [
    "mean_read_length",
    "median_read_length",
    "n50",
    "mean_quality",
    "median_quality",
    "total_bases"
]

for k in keys:
    report[k] = {
        "before": before.get(k),
        "after": after.get(k),
        "change": diff(before.get(k), after.get(k))
    }

# Quality thresholds
for q in ["Q10", "Q15", "Q20"]:
    key = f"reads_above_{q}"
    if key in before and key in after:
        report[key] = {
            "before": before[key]["percent"],
            "after": after[key]["percent"],
            "change": diff(before[key]["percent"], after[key]["percent"])
        }

# Save report
print(json.dumps(report, indent=2))