import pandas as pd
import numpy as np
from fancyimpute import SoftImpute
from sklearn.preprocessing import MinMaxScaler
import os

np.seterr(divide="ignore", invalid="ignore")
pd.options.mode.chained_assignment = None

def load_methylation_data(filepath):
    print("Loading raw methylation data...")
    raw = pd.read_csv(filepath, sep="\t", header=None)
    meth = raw.iloc[:, 0::2]
    depth = raw.iloc[:, 1::2]
    prop = pd.DataFrame(meth.values / depth.values)
    print(f"{prop.shape[0]} sites × {prop.shape[1]} samples loaded.")
    return meth, depth, prop

def perform_site_qc(prop, depth, tissue_info_path, ucla_n=95):
    print("Running site-level QC…")
    ucla_prop, uq_prop = prop.iloc[:, :ucla_n], prop.iloc[:, ucla_n:]
    ucla_depth, uq_depth = depth.iloc[:, :ucla_n], depth.iloc[:, ucla_n:]

    low_miss = set(ucla_prop.index[ucla_prop.isna().sum(1) < 5]) & \
                set(uq_prop.index[uq_prop.isna().sum(1) < 5])
    high_var = set(ucla_prop.index[ucla_prop.var(1) > ucla_prop.var(1).quantile(0.05)]) & \
               set(uq_prop.index[uq_prop.var(1) > uq_prop.var(1).quantile(0.05)])
    high_depth = set(ucla_depth.index[ucla_depth.median(1) > 10]) & \
                 set(uq_depth.index[uq_depth.median(1) > 10])

    tissues = pd.read_csv(tissue_info_path, sep="\t", header=None)
    tissues.columns = ["chr","start","end","name"]
    tissues["tissue"] = tissues["name"].str.split("-").str[1]
    non_placenta = set(tissues.index[tissues["tissue"] != "placenta"])

    keep = list(low_miss & high_var & high_depth & non_placenta)
    print(f"{len(keep)} CpG sites retained.")
    return keep

def impute_and_scale_cohort(prop, meta, group_name):
    samples = meta.loc[meta["group"] == group_name, "sample_num"]
    sub = prop[samples]

    imputer = SoftImpute()
    imputed = pd.DataFrame(imputer.fit_transform(sub), index=sub.index, columns=sub.columns)

    scaler = MinMaxScaler()
    scaled = pd.DataFrame(scaler.fit_transform(imputed.T).T, index=imputed.index, columns=imputed.columns)
    return imputed, scaled

def make_covariates(meta):
    from sklearn.preprocessing import OneHotEncoder, StandardScaler

    cov = meta[["sample_num", "Sex", "Ethnicity", "age", "concentration", "input_cfdna"]].copy()

    # One-hot encode categorical vars
    try:
        enc = OneHotEncoder(sparse_output=False, drop="first")  # sklearn >=1.2
    except TypeError:
        enc = OneHotEncoder(sparse=False, drop="first")         # sklearn <1.2

    cat = enc.fit_transform(cov[["Sex", "Ethnicity"]])

    try:
        cat_cols = enc.get_feature_names_out(["Sex", "Ethnicity"])  # sklearn >=1.0
    except AttributeError:
        cat_cols = enc.get_feature_names(["Sex", "Ethnicity"])      # sklearn <1.0

    cov_enc = pd.DataFrame(cat, columns=cat_cols, index=cov.index)

    # Scale numeric covariates
    cov_num = cov[["age", "concentration", "input_cfdna"]]
    scaler = StandardScaler()
    cov_num_scaled = pd.DataFrame(scaler.fit_transform(cov_num),
                                  columns=cov_num.columns, index=cov.index)

    cov_final = pd.concat([cov[["sample_num"]], cov_enc, cov_num_scaled], axis=1)
    return cov_final

def main():
    meth_path = "input_data.txt"
    meta_path = "metadata.csv"
    tissue_path = "sites.bed"

    outdir = "lasso_als_binary"
    os.makedirs(outdir, exist_ok=True)

    meth, depth, prop = load_methylation_data(meth_path)
    keep_sites = perform_site_qc(prop, depth, tissue_path)
    prop = prop.loc[keep_sites]

    meta = pd.read_csv(meta_path)
    prop.columns = meta["sample_num"]

    # Remove poor-quality samples
    missing_frac = prop.isna().sum() / len(prop)
    bad = list(missing_frac[missing_frac > 0.3].index)
    meta = meta[~meta["sample_num"].isin(bad)]
    prop = prop[[s for s in prop.columns if s not in bad]]

    meta = meta.dropna(subset=["age", "concentration", "input_cfdna", "Sex", "Ethnicity"])
    covs = make_covariates(meta)

    # cohort-specific imputation and scaling
    ucla_imp, ucla_scaled = impute_and_scale_cohort(prop, meta, "UCLA")
    uq_imp, uq_scaled     = impute_and_scale_cohort(prop, meta, "UQ")

    # attach covariates
    ucla_cov = covs[covs["sample_num"].isin(meta.loc[meta["group"]=="UCLA","sample_num"])]
    uq_cov   = covs[covs["sample_num"].isin(meta.loc[meta["group"]=="UQ","sample_num"])]

    ucla_scaled.to_csv(f"{outdir}/ucla_scaled.txt", sep="\t", index=False)
    uq_scaled.to_csv(f"{outdir}/uq_scaled.txt", sep="\t", index=False)
    ucla_cov.to_csv(f"{outdir}/ucla_covariates.txt", sep="\t", index=False)
    uq_cov.to_csv(f"{outdir}/uq_covariates.txt", sep="\t", index=False)

    meta[["sample_num","sample_type","group"]].to_csv(f"{outdir}/metadata_final.txt", sep="\t", index=False)
    
if __name__ == "__main__":
    main()
