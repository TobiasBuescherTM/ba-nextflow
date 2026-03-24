import pandas as pd
from open_nipals.nipalsPCA import NipalsPCA
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import os
from datetime import datetime
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("--input_data", required=True)
parser.add_argument("--file_name", required=True)
parser.add_argument("--log_file", required=True )
parser.add_argument("--log_output", required=True )
parser.add_argument("--output", required=True)
parser.add_argument("--pca_output", required=True)
args = parser.parse_args()

def qcomics (df): #missing values must be marked to use
    """
    -any missing values in an metabolite with more then 70% missing values is marked as a real Zero per group
    -and for values that are all missing in a biological group
    """
    
    bio_groups = set(df.loc["Group"].tolist())
    
    df_ = df.copy()
    drop_list = (["TECH_replica","BIO_replica","qc"])
    df_ = df_.drop(index = [index for index in drop_list if index in df_.index]  ) #drop row if avaiable
    df_ = df_.drop(columns = [col for col in ["IS"] if "IS" in df_.columns]) #drop column if avaiable 
    """
    if "qc" in df_.index:
        df_ = df_.drop(columns = df_.columns[df_.loc["qc"] == "TRUE"]) #drop qc samples 
    """

    treshhold = 0.7 #metabolites with 70% or more missing values are expected to be real zeros 
    real_zero_metabolites = []

    for group in bio_groups:

        group_columns = df_.columns[df_.loc["Group"] == group]
        df_group = df_[group_columns]

        for row in df_.index :

            if row == "Group":
                continue
            if row == "qc":
                continue

            counter = 0

            for col in df_group.columns:
                if df_.loc[row,col] == "-":
                    counter += 1 
    
            if counter / len(df_group.columns) >= treshhold:
                df_group_columns = df.columns[df.loc["Group"] == group]
                df.loc[row,df_group_columns] = df.loc[row,df_group_columns].replace("-", 0)
                real_zero_metabolites.append(row)

        

    
    return (df, real_zero_metabolites)
def NIPALS_PCA (df,path):
    """
    -PCA that skipes over missing values
    -to evaluate information gained /skewed with missing alue strategie based on
        -silluete score of group cluster

    """
    

    if "Group" not in df.index :
        raise ValueError("Groups not set ")
    
    if "qc" not in df.index :
        raise ValueError("QC Samples not set ")

    drop_list = (["TECH_replica"])
    df = df.drop(index = [index for index in drop_list if index in df.index]  ) #drop row if avaiable
    df = df.drop(columns = [col for col in ["IS"] if "IS" in df.columns]) #drop column if avaiable 
    
    Pca = NipalsPCA()
    metadata = df.loc[["Group","qc"]].T
    features_DF = df.drop(index=["Group","qc"]).T
    features = features_DF.to_numpy(dtype=float)
    
    pca = Pca.fit_transform(features)
    pca = pd.DataFrame(pca, columns=["PC1", "PC2"], index=features_DF.index)
    pca = pd.concat([pca, metadata], axis=1)

    qc = pca[pca["qc"] == "TRUE"][["PC1","PC2"]]
    centroid = qc.mean(axis=0)
    radius = np.max(np.linalg.norm(qc - centroid, axis=1))
    
    #print(pca)
    plt.figure(figsize=(10,7))
    sns.scatterplot(
        data=pca,
        x="PC1",
        y="PC2",
        hue="Group",       # color by group
        style="qc",        # shape by QC
        s=100
    )
    circle_before = Circle(centroid, radius, color="red", fill=False, linestyle="--", linewidth=2)
    plt.gca().add_patch(circle_before)
    plt.title("NIPALS-PCA")
    plt.xlabel(f"PC1")
    plt.ylabel(f"PC2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(path,f"nipals-pca.png"), dpi=300)
    plt.close()
def meancenter (df_):
    df = df_.copy()
    drop_list = (["TECH_replica","BIO_replica","qc","Group"])
    df = df.drop(index = [index for index in drop_list if index in df.index]  ) #drop row if avaiable
    df = df.drop(columns = [col for col in ["IS"] if "IS" in df.columns]) #drop column if avaiable 
    
    for row in df.index :
        row_data = df.loc[row].tolist()
        row_mean = np.nanmean([float(x) for x in row_data if x != "-"])
        new_row = []
        for value in row_data:
            if pd.isna(value) or value == "-":
                new_row.append("-")
            else:
                new_row.append(float(value) - row_mean)

        df_.loc[row] = new_row
    
    return df_   
def zscale (df):
    df_ = df.copy()
    drop_list = (["Group","TECH_replica","BIO_replica","qc"])
    df_ = df_.drop(index = [index for index in drop_list if index in df_.index]  ) #drop row if avaiable
    df_ = df_.drop(columns = [col for col in ["IS"] if "IS" in df_.columns]) #drop column if avaiable

    for row in df_.index:

        values = df_.loc[row].tolist()
        values_ = [float(x) for x in values if x != "-"]

        values_mean = np.mean(values_) 
        values_sd = np.std(values_)

        zscores = []
        for i in range(len(values)):
           
           
           if values[i] == "-":
               zscores.append("-")
               continue

           zscores.append((float(values[i]) - values_mean) / values_sd)

        df.loc[row] = zscores 

        
    return df
def NIPALS_PCA_multible (df,path):
    """
    -PCA that skipes over missing values
    -to evaluate information gained /skewed with missing alue strategie based on
        -silluete score of group cluster

    """
    

    if "Group" not in df.index :
        raise ValueError("Groups not set ")
    
    if "qc" not in df.index :
        raise ValueError("QC Samples not set ")

  
    df = df.drop(columns = [col for col in ["IS"] if "IS" in df.columns]) #drop column if avaiable 
    
    Pca = NipalsPCA()
    metadata = df.loc[["Group","qc","TECH_replica"]].T
    features_DF = df.drop(index=["Group","qc","TECH_replica"]).T
    features = features_DF.to_numpy(dtype=float)
    
    pca = Pca.fit_transform(features)
    pca = pd.DataFrame(pca, columns=["PC1", "PC2"], index=features_DF.index)
    pca = pd.concat([pca, metadata], axis=1)

    print(pca)



    qc = pca[pca["qc"] == "TRUE"][["PC1","PC2"]]
    centroid = qc.mean(axis=0)
    radius = np.max(np.linalg.norm(qc - centroid, axis=1))
    unique_tech_replicas = sorted(set(x for x in pca["TECH_replica"].tolist() if x != "FALSE"))
    n_replicates = len(unique_tech_replicas)
    
    n_plots = len(unique_tech_replicas) + 1  # +1 for main PCA plot
    n_cols = min(3, n_replicates)            # max 3 per row
    n_rows_reps = math.ceil(n_replicates / n_cols)

    # create figure with GridSpec
    fig = plt.figure(figsize=(6*n_cols, 5*(1 + n_rows_reps)))
    gs = fig.add_gridspec(1 + n_rows_reps, n_cols, height_ratios=[2] + [1]*n_rows_reps)

    # main PCA plot spans top row
    ax = fig.add_subplot(gs[0, :])















    
    
    sns.scatterplot(
        data=pca,
        x="PC1",
        y="PC2",
        hue="Group",       # color by group
        style="qc",        # shape by QC
        s=100,
        ax = ax
    )

    
    circle_before = Circle(centroid, radius, color="red", fill=False, linestyle="--", linewidth=2)
    ax.add_patch(circle_before)
    ax.set_title("NIPALS-PCA")
    ax.set_xlabel(f"PC1")
    ax.set_ylabel(f"PC2")
    ax.set_aspect('equal', adjustable='datalim')
    #plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    #plt.tight_layout()
    #plt.savefig(path, dpi=300)
    #plt.show(block = True)
    #plt.close()






    axes_reps = []

    for idx, replica in enumerate(unique_tech_replicas):
        row = 1 + idx // n_cols       # row in gridspec
        col = idx % n_cols            # column in gridspec
        ax = fig.add_subplot(gs[row, col])
        replicates = pca[pca["TECH_replica"] == replica][["PC1","PC2"]]
        centroid = replicates.mean(axis=0)
        radius = np.max(np.linalg.norm(replicates - centroid, axis=1))

        replicates = pca[pca["TECH_replica"] == replica][["PC1","PC2","Group"]]
        


       
        sns.scatterplot(
            data=pca,
            x="PC1",
            y="PC2",
            hue="Group", 
            alpha = 0.3,
            s=40,
            ax=ax
            
        )

        group_name = replicates["Group"].iloc[0]
        color = ax.get_legend().get_lines()[list(pca["Group"].unique()).index(group_name)].get_color()

        sns.scatterplot(
            data = replicates,
            x="PC1",
            y="PC2",
            facecolor="none",     
            edgecolor=color,      
            linewidth=2,
            s=120,
            ax = ax,
            legend = False
        )

        title_string = "replicate samples "
        for j, row in replicates.iterrows():
            ax.text(
                row["PC1"] + 0.02,   
                row["PC2"] + 0.02,
                str(row.name),        
                fontsize=12,
                color="black"
            )
            title_string =  title_string + str(row.name) + " "
        


        circle = Circle(centroid, radius, color=color, fill=False, linestyle="--", linewidth=2)
        ax.add_patch(circle)
        
        ax.set_title(title_string)
        ax.set_xlabel(f"PC1")
        ax.set_ylabel(f"PC2")
        ax.set_aspect("equal", adjustable="datalim")
        
       
    plt.tight_layout()  
    #plt.show(block = True)
    plt.savefig(os.path.join(path,f"nipals_pca_with_replicates.png"), dpi=300)
    plt.close()

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

with open(args.log_file, "r") as src, open(args.log_output, "w") as dst:
    dst.write(src.read())

with open(args.log_output, "a") as file:
    file.write(f"qcomicsStart: {timestamp}\n")

os.makedirs(args.pca_output, exist_ok=True)


df_in = pd.read_csv(args.input_data, sep = ",", index_col=0)
df_out,metabolites_effected = qcomics(df_in)
df_out.to_csv(args.output, sep=",")

df_centered = meancenter(df_out)
df_scaled = zscale(df_centered)
df_replaced = df_scaled.replace("-",np.nan)
NIPALS_PCA(df_replaced,args.pca_output)
NIPALS_PCA_multible(df_replaced,args.pca_output)

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

with open(args.log_output, "a") as file:
    file.write(f"metabolites_effected: {metabolites_effected}\n")
    file.write(f"qcomicsClose: {timestamp}\n")