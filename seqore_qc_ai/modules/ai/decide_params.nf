process DECIDE_FILTER_PARAMS {

    input:
    path metrics_json

    output:
    val params

    script:
    """
    python3 - <<EOF
import json

with open("${metrics_json}") as f:
    m = json.load(f)

q15 = m.get("reads_above_Q15", {}).get("percent", 0)

# Simple rule logic
if q15 < 40:
    params = {"min_q": 12, "min_len": 2000}
else:
    params = {"min_q": 10, "min_len": 1000}

print(params)
EOF
    """
}