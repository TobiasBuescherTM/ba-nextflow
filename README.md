
## Installation

First, clone the repository:

```bash
git clone https://github.com/TobiasBuescherTM/ba-nextflow.git
```

There are two ways to run the pipeline.

### Option 1: All-in-one Docker container

Everything runs inside a single container.

1. Pull the Docker image:

   ```bash
   docker pull tobiasbueschertm/ba-nf
   ```

2. Run the container and mount the cloned repository:

   ```bash
   docker run -it --rm \
     -v "$(pwd)":/workspace \
     tobiasbueschertm/ba-nf
   ```

3. Change into the mounted directory inside the container:

   ```bash
   cd /workspace
   ```

### Option 2: Multiple Docker containers (Linux required)

Each module runs in its own container, orchestrated by Nextflow.

1. Get the Docker images for all modules defined in the repository, either by pulling them from `https://hub.docker.com/u/tobiasbueschertm`:

   ```bash
   docker pull tobiasbueschertm/nf-missing_percent
   docker pull tobiasbueschertm/nf-cv_icc
   docker pull tobiasbueschertm/nf-imputation
   docker pull tobiasbueschertm/nf-innit
   docker pull tobiasbueschertm/nf-qcomics
   docker pull tobiasbueschertm/nf-replicate_zero
   docker pull tobiasbueschertm/nf-vsn
   ```

   or by building them from the Dockerfiles inside the repository:

   ```bash
   docker build -t tobiasbueschertm/nf-missing_percent -f containers/Dockerfile.80% .
   docker build -t tobiasbueschertm/nf-cv_icc -f containers/Dockerfile.cv_icc .
   docker build -t tobiasbueschertm/nf-imputation -f containers/Dockerfile.imputation .
   docker build -t tobiasbueschertm/nf-innit -f containers/Dockerfile.innit .
   docker build -t tobiasbueschertm/nf-qcomics -f containers/Dockerfile.qcomics .
   docker build -t tobiasbueschertm/nf-replicate_zero -f containers/Dockerfile.replicate_zero .
   docker build -t tobiasbueschertm/nf-vsn -f containers/Dockerfile.vsn .
   ```

2. Download and install Nextflow from `https://www.nextflow.io`, for example:

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. Make sure Nextflow is on your `PATH`:

   ```bash
   export PATH="$PATH:/path/to/nextflow/directory"
   ```

   Add this line to your shell config (e.g. `~/.bashrc`) to make it permanent.

4. Set the syntax parser variable, either for the current session:

   ```bash
   export NXF_SYNTAX_PARSER=v2
   ```

   or permanently by adding the same line to your shell config (e.g. `~/.bashrc` or `~/.zshrc`).

5. In `nextflow.config`, enable Docker:

   ```groovy
   docker.enabled = true
   ```

## Running an experiment

1. **Format your data.** Use `csv_formater.py` to load your MS data table, set the data groups (biological groups, technical variations and QC samples), and export the formatted table.

2. **Add the data.** Put the formatted table into the `input_data` folder.

3. **Run the pipeline.** From the root directory of the repository, execute:

   ```bash
   nextflow run main.nf -name "Experiment name"
   ```

4. **Collect the results.** Once all steps have completed, you will find all results in the `results/outputs` folder under the experiment name you set.

## Abstract

The objective of this work was to develop an automated preprocessing pipeline specifically tailored to a given MS metabolomics dataset generated from a case-control study comparing patient samples with a specific disease against a control pool of “healthy” individuals. Although several integrated pipelines and toolkits, such as Rodin, TidyMS, the nPYc Toolbox, SECIMTools, and pyOpeMS, are freely available, the given dataset contains specific characteristics, including pooled QC samples and technical replicates, that are not fully addressed in these existing solutions. These features can potentially enhance preprocessing performance by enabling more informed assumptions about the underlying biological truth. 
Its three complementary strategies are to classify missing values and impute them with MAI guided by technical-replicate behavior, apply CV-ICC and missing-abundance quality filtering to distinguish reliable from questionable signals, and normalize the data using a variance-stabilizing approach (VSN). Technical replicates help distinguish true blanks from values requiring imputation, while pooled QCs support more accurate assessment of metabolite data quality across the run. 
Proof-of-concept evaluation used two artificial datasets derived from the real dataset and engineered to reflect high versus reduced biological variation, PCA cluster quality and machine-learning classification performance were compared with and without the new missing-value strategy. Cluster separation improved clearly in both datasets, whereas classification performance increased modestly for the uncut dataset and declined more noticeably for the cut dataset. The workflow is considered validated as a concept, and the classification decline after filtering is attributed to the stricter 80% rule rather than to a methodological flaw. Because boolean missing-value classifications benefit downstream methods only when those methods can use them, a zero-missing strategy is recommended for methods that cannot. 
Limitations include the potential for imputation bias, such as imputing metabolites associated with a specific sex across groups, use of the original MAI method because the newer PSO/XGBoost variant lacks open-source code, reduced classification power under the 80% rule, for which projecting real-value distributions onto placeholder NAs is proposed as an alternative, variable QRILC performance, necessitating broader comparisons, and the lack of benchmarking of the Nextflow implementation. The results remain a proof of concept and are not directly comparable with other pipelines. The data generator could not model technical-replicate variation and assumed nearly identical, normal-like distributions across groups, consequently, some methods were not fully evaluated, and privacy constraints prevented use of the real dataset.
