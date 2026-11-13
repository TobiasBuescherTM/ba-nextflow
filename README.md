# Project Name

Short description of what this pipeline does.

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

2. Download and install Nextflow from `<NEXTFLOW_LINK>`, for example:

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

2. **Add the data.** Put the formatted table into the `input` folder.

3. **Run the pipeline.** From the root directory of the repository, execute:

   ```bash
   nextflow run main.nf -name "Experiment name"
   ```

4. **Collect the results.** Once all steps have completed, you will find all results in the `output` folder under the experiment name you set.
