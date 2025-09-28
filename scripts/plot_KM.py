import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
import re


def KM_plot_missing(df, col, col_event, col_time):
    kmf = KaplanMeierFitter()

    mask_missing = df[col].isna()
    time = df[col_time]
    event = df[col_event]

    plt.figure(figsize=(7,5))

    kmf.fit(time[~mask_missing], event[~mask_missing], label = col)
    kmf.plot_survival_function()

    kmf.fit(time[mask_missing], event[mask_missing], label = col)
    kmf.plot_survival_function()

    plt.title(f"Kaplan Meier curves - {col}")
    
    
    
def detect_monosomie7(karyotype):
    """
    Return 1 if monosomie 7 or del(7q), else 0.
    """
 
    if re.search(r"(?<!\d)-7(?!\d)", karyotype):
        return 1
    

    if re.search(r"del\(7q\d+\)", karyotype):
        return 1
    
    return 0