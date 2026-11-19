# qfit_pipeline configuration (bash syntax; sourced by every script).
#
# A run copies this file to <run_dir>/run_info/config.sh when it is created and
# always uses that copy afterwards, so editing this file never changes a run
# that already exists. Use a new run directory to try different settings.

#______________________________ SOFTWARE ______________________________
# Apptainer image with Phenix and the fixed qFit (see container/qfit.def).
# Leave empty to run with Phenix and qFit from the current environment
# instead (e.g. a workstation where `source phenix_env.sh` and
# `conda activate qfit` have been run).
CONTAINER="/path/to/qfit_pipeline.sif"

# Host paths visible inside the container. Must include this repository, the
# run directory and INPUT_DIR. Apptainer also binds $HOME and the working
# directory by default.
APPTAINER_BINDS="/path/to/data"

# Optional: abort if the qFit in the container is not this commit (the
# "+g<hash>" part of its version, or /opt/qfit-3.0.commit in the image).
QFIT_EXPECTED_COMMIT=""

#______________________________ INPUTS ______________________________
# Empty: download each structure with phenix.fetch_pdb.
# Otherwise a directory laid out as INPUT_DIR/<ID>/<id>.pdb and <id>.mtz
# (optional: <id>-sf.cif, <id>.fa), e.g. the test_data folder.
INPUT_DIR=""

#______________________________ PRE-QFIT REFINEMENT ______________________________
PRE_MACRO_CYCLES=8
# refinement.target_weights.optimize_xyz_weight / optimize_adp_weight.
# Pending the weight-optimization test; false matches the original pipeline.
OPTIMIZE_WEIGHTS=false

#______________________________ QFIT ______________________________
# qfit_protein --scale-rmask: scales the radius of the atom mask used to fit
# the map scale. With the default (1.0) the fit sees almost no solvent, which
# leaves the solvent region of the scaled map too high when model-map agreement
# is modest; qFit then adds spurious alternate conformations. 2.0 brought the
# solvent level to ~0.33 in sampling tests on 28 structures, removing most
# spurious alternates, but has not yet been validated by refinement. 1.0 is
# qFit's default.
QFIT_SCALE_RMASK=1.0

#______________________________ RESOURCES ______________________________
# Cores, memory and time limit per stage. On SLURM these become
# --cpus-per-task, --mem and --time; locally only the cores are used.
PRE_CPUS=1
PRE_MEM=4G
PRE_TIME=10:00:00

QFIT_CPUS=20
QFIT_MEM=18G
QFIT_TIME=10:00:00

POST_CPUS=1
POST_MEM=4G
POST_TIME=100:00:00

#______________________________ SLURM ______________________________
SLURM_ACCOUNT=""       # --account, if your cluster needs one
SLURM_PARTITION=""     # --partition, if your cluster needs one
SLURM_EXTRA_ARGS=""    # anything else for sbatch, e.g. "--qos=long"

# Maximum number of array tasks of each stage running at once.
LIMIT_PRE=100
LIMIT_QFIT=40
LIMIT_POST=100
