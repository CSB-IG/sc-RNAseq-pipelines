import scanpy 
import celltypist 
from celltypist import models

sc.settings.set_figure_params(dpi=50, facecolor="white")


path_results = "your/results/"
if not os.path.exists(path_results):
    os.makedirs(path_results, exist_ok=True)

path_file = "./preprocess_FILE.h5ad"

adata = sc.read_h5ad(path_file)

models.download_models(force_update = True)

models.models_description()

model = models.Model.load(model = 'YOUR_MODEL.pkl')

predictions = celltypist.annotate(adata, model = model, majority_voting = True)

predictions.predicted_labels

adata = predictions.to_adata()

sc.tl.umap(adata)

sc.pl.umap(adata, color = ['cell_type', 'predicted_labels', 'majority_voting'], legend_loc = 'on data')

celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting')

adata.write_h5ad("./auto_annot.py")
