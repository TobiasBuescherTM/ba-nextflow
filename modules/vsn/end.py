import numpy as np
import math
import os
from datetime import datetime
import argparse
import pandas as pd
from open_nipals.nipalsPCA import NipalsPCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle
from sklearn.decomposition import PCA

parser = argparse.ArgumentParser()

parser.add_argument("--file_name", required=True)
parser.add_argument("--log_output", required=True )
parser.add_argument("--output", required=True)
parser.add_argument("--pca_output", required=True)
args = parser.parse_args()

os.makedirs(args.pca_output, exist_ok=True)


def PCA_func (df,path):
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
    
    Pca = PCA()
    metadata = df.loc[["Group","qc"]].T
    features_DF = df.drop(index=["Group","qc"]).T
    features = features_DF.to_numpy(dtype=float)
    
    pca_fit = Pca.fit_transform(features)
    pca_obj = Pca.fit(features)

    explained_var = pca_obj.explained_variance_ratio_
    components = np.arange(1, len(explained_var)+1)
    top2_idx = np.argsort(explained_var)[-2:][::-1]
    pca_top2 = pca_fit[:, top2_idx]


    pca = pd.DataFrame(pca_top2, columns=["PC1", "PC2"], index=features_DF.index)
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
    plt.title("PCA")
    plt.xlabel(f"PC1")
    plt.ylabel(f"PC2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(path,f"pca.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(10,6))
    plt.bar(components, explained_var*100, alpha=0.7, color="skyblue", label="Individual")
    plt.plot(components, np.cumsum(explained_var)*100, marker="o", color="red", label="Cumulative")
    plt.xlabel("Principal Component")
    plt.ylabel("Explained Variance (%)")
    plt.title("Explained Variance by PCA Components with Normalization")
    plt.xticks(components)
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(path,f"pca_explained_var.png"), dpi=300)
    plt.close()
def PCA_multible (df,path):
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
    
    Pca = PCA()
    metadata = df.loc[["Group","qc","TECH_replica"]].T
    features_DF = df.drop(index=["Group","qc","TECH_replica"]).T
    features = features_DF.to_numpy(dtype=float)
    
    pca_obj = Pca.fit(features)
    pca = Pca.fit_transform(features)
    explained_var = pca_obj.explained_variance_ratio_
    components = np.arange(1, len(explained_var)+1)
    top2_idx = np.argsort(explained_var)[-2:][::-1]
    pca_top2 = pca[:, top2_idx]


    pca = pd.DataFrame(pca_top2, columns=["PC1", "PC2"], index=features_DF.index)
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
    ax.set_title("PCA")
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
    plt.savefig(os.path.join(path,f"pca_with_replicates.png"), dpi=300)
    plt.close()


df_out = pd.read_csv(args.output, sep = ",", index_col=0)

PCA_func(df_out,args.pca_output)
PCA_multible(df_out,args.pca_output)

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

with open(args.log_output, "a") as file:
    file.write(f"vsnClose: {timestamp}\n")