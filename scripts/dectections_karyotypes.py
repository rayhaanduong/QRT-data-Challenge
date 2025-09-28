import re
import re
from scipy.stats import chi2_contingency, fisher_exact
import pandas as pd
from collections import Counter


def detect_monosomie7(karyotype):
    """
    Return 1 if monosomie 7 or del(7q), else 0.
    """
 
    if re.search(r"(?<!\d)-7(?!\d)", karyotype):
        return 1
    

    if re.search(r"del\(7q\d+\)", karyotype):
        return 1
    
    return 0



def detect_anomaly(karyotype, anomaly_pattern):
    anomaly_pattern = re.escape(anomaly_pattern)  
    return int(bool(re.search(anomaly_pattern, karyotype, re.IGNORECASE)))

def cytogenetic_analysis(df, anomalies, target_col="TARGET"):
    """
    
    Args: 
        df : pd.DataFrame : dataframe with CYTOGENETICS and TARGET.
        target_col : str
        anomalies : List : List of anomalies
    
    Return Dataframe with one-hot encoding + p-values.
    """
    results = []
    
    for anomaly in anomalies:
        df[anomaly] = df["CYTOGENETICS"].apply(lambda x: detect_anomaly(str(x), anomaly))
        
        contingency = pd.crosstab(df[anomaly], df[target_col])
        
        
        if contingency.shape == (2, 2) and contingency.values.min() < 5:
            _, p = fisher_exact(contingency)
        else:
            _, p, _, _ = chi2_contingency(contingency)
        
        results.append({"anomaly": anomaly, "p_value": p, "present_cases": df[anomaly].sum()})
    
    results_df = pd.DataFrame(results).sort_values("p_value")
    return df, results_df




def find_anomalies(df, k = 50):
    cyto_series = df["CYTOGENETICS"].dropna().astype(str)


    pattern = r"(\-?\+?\d+|del\(\d+q?\)|del\(\d+p?\)|t\(\d+;\d+\)|inv\(\d+\)|complex karyotype)"

    all_matches = []
    for entry in cyto_series:
        matches = re.findall(pattern, entry, re.IGNORECASE)
        all_matches.extend(matches)


    freq_anomalies = Counter(all_matches)


    freq_df = pd.DataFrame(freq_anomalies.items(), columns=["anomaly", "count"]).sort_values(by="count", ascending=False)
    
    return freq_df[freq_df["count"] >= k]["anomaly"].to_list()