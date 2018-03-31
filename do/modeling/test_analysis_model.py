
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.analysis.run import AnalysisRuns
from pts.core.basics.log import setup_log
from pts.core.plot.sed import plot_sed

setup_log("DEBUG")

modeling_path = verify_modeling_cwd()

runs = AnalysisRuns(modeling_path)

run = runs.load("new_lowres2")

model = run.model

sed = model.intrinsic_sed_sfr

# Plot
plot_sed(sed, "SFR")
