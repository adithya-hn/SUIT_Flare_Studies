import pandas as pd

df_peaks = pd.read_csv("suit_diff_peaks.csv", parse_dates=["date_time"])
df_flux  = pd.read_csv("v2_Diff_img_data_NB04.csv", parse_dates=["date_time"])
df_flux_pos = df_flux[df_flux["diff_count"] > 0]

df_common = pd.merge_asof(
    df_peaks.sort_values("date_time"),
    df_flux_pos.sort_values("date_time"),
    on="date_time",
    tolerance=pd.Timedelta("30s"),
    direction="nearest"
)

df_common.to_csv("c10_common_peak_list_flux_gt0.csv", index=False)
print("Saved common_peak_list_flux_gt0.csv")