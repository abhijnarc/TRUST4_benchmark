python3 - <<'PY'
import pandas as pd

df = pd.read_csv("PRJNA492301_runinfo.csv")

print("Number of runs:", len(df))
print("\nColumns:")
print(df.columns.tolist())

print("\nRuns:")
cols = [c for c in ["Run", "SampleName", "LibraryName",
                    "LibraryLayout", "spots", "bases", "Platform"]
        if c in df.columns]

print(df[cols].to_string(index=False))
PY