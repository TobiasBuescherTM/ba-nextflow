import pandas as pd
from open_nipals.nipalsPCA import NipalsPCA
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import os
from datetime import datetime
import argparse
import pingouin as pg
from scipy.stats import pearsonr
import os
import math


parser = argparse.ArgumentParser()
parser.add_argument("--input_data", required=True)
parser.add_argument("--file_name", required=True)
parser.add_argument("--cv_cutoff",type=float, required=True)
parser.add_argument("--icc_cutoff",type=float, required=True)
parser.add_argument("--log_file", required=True )
parser.add_argument("--log_output", required=True )
parser.add_argument("--output", required=True)
parser.add_argument("--pca_output", required=True)
parser.add_argument("--correlations_path", required=True)
args = parser.parse_args()



        

def cv_filter (df,treshhold):

    """
    -detect all metabolites with an CV over a specified values
    """

    if "qc" not in df.index:
        raise ValueError("QC-samples not set ")
    
    
    df_ = df.copy()
    drop_list = (["Group","TECH_replica","BIO_replica"])
    df_ = df_.drop(index = [index for index in drop_list if index in df_.index]  ) #drop row if avaiable
    df_ = df_.drop(columns = [col for col in ["IS"] if "IS" in df_.columns]) #drop column if avaiable 

   
    filtered_metabolites = []
    data = []
    treshhold = treshhold #treshold not clear 0.2-0.3


    df_ = df_.drop(columns = df_.columns[df_.loc["qc"] == "FALSE"]) #drop all no qc-samples 
    #print (df_)

    for row in df_.index:
        if row == "qc":
            continue

        data = df_.loc[row].tolist()
        data = [float(value) for value in data if value != "-" ] #filter out all missing values 
        data = [value for value in data if value != 0 ]
        #print(data)
        if len(data) < 3 :
            print("metabolite has to much missing values (metabolite id filtered out [CV-filter])")
            filtered_metabolites.append(row)
            continue

        cv = np.std(data, ddof=1) / np.mean(data) #using 1 degree of freedom / sample standart deviation

        if np.isnan(cv):
            print("NaN detcted in Cv filter")
            print(data)
       

        if cv > treshhold:
            filtered_metabolites.append(row)
            



    


   



    return (filtered_metabolites)
def ICC_filter (df,metabolite_list,treshhold): #ICC with one way random effects model #using pingouin to calculate 

    """
    -filter the detected metabolite from the ICC_filter 
    -if ICC is below a specified value the metabolite is deleted
    """

    
    filtered_metabolites = []

    

    def metabolite_icc(flat_df):

        icc_table = pg.intraclass_corr(data=flat_df, targets="sample",raters="replicate",ratings="intensity")

        icc = icc_table.loc[icc_table["Type"] == "ICC1","ICC"].item()    #one way randome effect

        return icc


    if "TECH_replica" not in df.index:
        raise ValueError("technical replicas not set ")
    
    
    df_ = df.copy()
    drop_list = (["Group","BIO_replica","qc"]) 
    df_ = df_.drop(index = [index for index in drop_list if index in df_.index]  ) #drop row if avaiable
    df_ = df_.drop(columns = [col for col in ["IS"] if "IS" in df_.columns]) #drop column if avaiable 
    df_ = df_.drop(columns = df_.columns[df_.loc["TECH_replica"] == "FALSE"]) #drop all samples with no replicas

    
    for row in metabolite_list: #df_.index for testing

        if row == "TECH_replica":
            continue

        flat_df = pd.DataFrame(columns=["sample","replicate","intensity"])
        sample = []
        replicate = []
        intensity = []
        test_list = []
        saved_metabolites = []

        values = df_.loc["TECH_replica"].tolist()
        #print(df_)

        for value in set(values):

            test_list.clear()
            test_list = df_.loc[row, df_.columns[df_.loc["TECH_replica"] == value]].tolist()
            if "-" in test_list: #if replicate contains missing value 
                continue

            

            intensity += np.array(test_list,dtype=float).tolist()

            intensity
            sample += ["s" + str(value)] * values.count(value)

            for i  in range(1,values.count(value)+1): #skip zero
                replicate.append(i)
           
            
        #print(sample,replicate,intensity)
        flat_df["sample"] = sample
        flat_df["replicate"] = replicate
        flat_df["intensity"] = intensity

        icc = metabolite_icc(flat_df)
        if icc < treshhold :
            filtered_metabolites.append(row)
        else:
            saved_metabolites.append(row)

    df = df.drop (index = [index for index in filtered_metabolites])
    return (df,filtered_metabolites,saved_metabolites)
def CV_Intensiety_Correlation (df,ICC_METABOLITE_LIST,CV_TRESHHOLD):
    
    
    """
    -corrolate the mean intensiety vs the mean CV between technical replicas for every biological group
    -Evaluate CV and ICC Filter system 
    """

    if "TECH_replica" not in df.index:
        raise ValueError("Tech-Replica not set ")
    
    if "Group" not in df.index:
        raise ValueError("Group not set ") 
    

    os.makedirs(args.correlations_path, exist_ok=True)
    open(os.path.join(args.correlations_path,f"/imHereToAvoidAnErrorIgnoreMe.txt"),"w").close()
    
    df_reformat = df.copy()
    drop_list = (["qc","BIO_replica"])
    df_reformat = df_reformat.drop(index = [index for index in drop_list if index in df_reformat.index]  ) #drop row if avaiable
    df_reformat = df_reformat.drop(columns = [col for col in ["IS"] if "IS" in df_reformat.columns]) #drop column if avaiable 
    df_reformat = df_reformat.loc[:, df_reformat.loc["TECH_replica"] != "FALSE"] #exclude all rows with no tech-replicas

    group_data = {}
    correlation_data = {}
    bio_group_without_correlation = []



    for biogroup in set(df_reformat.loc["Group"].tolist()):

        try:
            excludet_mbt = []
            
            df_ = df_reformat
            df_ = df_[df_.columns[df_.loc["Group"] == biogroup]]
            #print(df_)

            for row in df_.index:
                if row == "TECH_replica":
                    continue

                if row == "Group":
                    continue
                
                replica_cvs = []
                for replica_set in set(df_.loc["TECH_replica"].tolist()):

                    values = df_.loc[row, df_.columns[df_.loc["TECH_replica"] == replica_set]].tolist()
                    values = [float(x) for x in values if  x != "-" ] #use only real values no zero ore missing values
                    values = [x for x in values if x != 0]

                    if len(values) < 2:
                        continue

                    cv = np.std(values, ddof=1) / np.mean(values)#using 1 degree of freedom / sample standart deviation)

                    if np.isnan(cv):
                        print("NAN value detected")
                        print(values)

                    replica_cvs.append(cv) #using 1 degree of freedom / sample standart deviation)

                all_values = df_.loc[row].tolist()

                all_values = [float(x) for x in all_values if  x != "-" ]
                all_values = [x for x in all_values if x != 0]

                if len(all_values) < 2:
                    excludet_mbt.append(row)
                    continue
                
                
                
                mean = np.mean(all_values)
                cv_means = np.mean(replica_cvs) 
                
                if len(replica_cvs) < 2:
                    excludet_mbt.append(row)
                    continue

                correlation_data[row] = [cv_means,mean]

                


            marked_mean_data = []
            marked_cv_data = []
            correlation_mean_data = []
            correlation_cv_data = []

            #print(correlation_data)

            for row,data_list in correlation_data.items():
                if row in ICC_METABOLITE_LIST:
                    marked_mean_data.append(data_list[1])
                    marked_cv_data.append(data_list[0])
                    continue

                correlation_mean_data.append(data_list[1])
                correlation_cv_data.append(data_list[0])

            #print(correlation_mean_data)
            #print(correlation_cv_data)

            corr, p_value = pearsonr(correlation_mean_data, correlation_cv_data)

            group_data[biogroup] = [corr, p_value]

            #scatterplot with best fit linear line
            correlation_mean_data = np.array(correlation_mean_data)
            correlation_cv_data = np.array(correlation_cv_data)

            

            m, b = np.polyfit(correlation_mean_data, correlation_cv_data, 1)  # linear fit

            fig = plt.figure(figsize=(8, 6))
            plt.scatter(correlation_mean_data, correlation_cv_data)
            plt.plot(correlation_mean_data, m*correlation_mean_data + b)
            plt.scatter(marked_mean_data,marked_cv_data,s=120, marker='x')#highlight filtered metabolites 
            #label filtered metabolites
            for xi, yi, label in zip(marked_mean_data, marked_cv_data, ICC_METABOLITE_LIST):
                plt.annotate(label,(xi, yi),textcoords="offset points",xytext=(5, 5),fontsize=5)
            plt.axhline(CV_TRESHHOLD, linestyle="--", linewidth=2, color="red", label="Filter - Threshold") #filter line
            plt.xlabel("mean")
            plt.ylabel("CV")
            plt.title("Correlation of mean vs Coefficient  \n of Variation with Regression Line \n Biogroup: " + str(biogroup))
            #plt.xlim(0, 1e6)
            plt.legend()
            excluded_text = "Excluded metabolites:\n" + ", ".join(excludet_mbt)
            plt.subplots_adjust(bottom=0.2)
            fig.text(0.5, 0.02, excluded_text, ha="center", va="bottom", fontsize=9)
            plt.savefig(os.path.join(args.correlations_path,f"cv_mean_correlation_biogroup_{biogroup}.png"), dpi=300)
            plt.close()
        except:
            bio_group_without_correlation.append(biogroup)


    return (group_data,bio_group_without_correlation)
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
    file.write(f"cv_iccStart: {timestamp}\n")

os.makedirs(args.pca_output, exist_ok=True)


df_in = pd.read_csv(args.input_data, sep = ",", index_col=0)
cv_metabolites = cv_filter(df_in,args.cv_cutoff)
df_out,cut_metabolites,metabolites_above_icc = ICC_filter(df_in,cv_metabolites,args.icc_cutoff) 
group_correlations,bio_group_without_correlation  = CV_Intensiety_Correlation(df_in,cut_metabolites,args.cv_cutoff)







df_out.to_csv(args.output, sep=",")

df_centered = meancenter(df_out)
df_scaled = zscale(df_centered)
df_replaced = df_scaled.replace("-",np.nan)
NIPALS_PCA(df_replaced,args.pca_output)
NIPALS_PCA_multible(df_replaced,args.pca_output)

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

with open(args.log_output, "a") as file:
    file.write(f"metabolites above cv cutoff: {cv_metabolites}\n")
    file.write(f"metabolites above icc cutoff: {metabolites_above_icc}\n")
    file.write(f"metabolites filtered out: {cut_metabolites}\n")
    file.write(f"cv vs mean_intensity correlation per bio group: {group_correlations}\n")
    file.write(f"bio group without correlation data (most likely to many missing values): {bio_group_without_correlation}\n")
    file.write(f"cv_iccClose: {timestamp}\n")

