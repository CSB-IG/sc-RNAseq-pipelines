#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('cd', '/home/mdiaz/hackaton-csbig/ranger/piloto/clean_adata/')


# In[1]:


import numpy as np
import scanpy as sc
import seaborn as sns
import pandas as pd 
import scvi
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
import os
from scipy.stats import median_abs_deviation as mad
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=False,
)


# In[3]:


os.listdir ("../clean_adata/")


# # QC preprocess for leukimia dataset 

# In[15]:


adatas = [x for x in os.listdir('../clean_adata/') if x.endswith('filtered.h5')]
adatas


# In[16]:


adatas


# In[17]:


def load_it(adata):
    samp = adata.split('_')[0]
    dx = adata.split('_')[1]
    adata = sc.read_10x_h5('../clean_adata/' + adata)
    adata.obs['Patient'] = samp
    adata.obs['DX'] = dx
    adata.obs['Sample'] = adata.obs['Patient'] + '_' + adata.obs['DX']
    adata.obs.index = adata.obs.index + '-' + samp + '_' + dx
    return adata


# In[19]:


adatas = [load_it(ad) for ad in adatas]

for adata in adatas:
    adata.var_names_make_unique()
    print(f"Unique var names for {adata.obs['Sample'][0]}: {adata.var_names.is_unique}")

adatas


# In[25]:


def qc(adata):
    
    sc.pp.filter_cells(adata, min_genes = 200) 
    
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars="mt", inplace=True, percent_top=[20])
    sc.pp.normalize_total(adata, target_sum=10**4)
    sc.pp.log1p(adata)

    remove = ['total_counts_mt', 'log1p_total_counts_mt']
    
    adata.obs = adata.obs[[x for x in adata.obs.columns if x not in remove]]
    
    return adata


# In[26]:


adatas = [qc(ad) for ad in adatas]


# In[27]:


df = pd.concat(x.obs for x in adatas)


# In[28]:


df = df.sort_values('Sample')


# In[29]:


df


# In[30]:


#value = "pct_counts_mt"
value = "n_genes"
#value = 'pct_counts_in_top_20_genes'
#value = "log1p_total_counts"

sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

g = sns.FacetGrid(df, row="Sample", hue="Sample", aspect=10, height=1, palette="tab20")

g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

g.map(plt.axhline, y=0, lw=2, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)

g.map(label, value)

g.figure.subplots_adjust(hspace=-0.2)

g.set_titles("")
g.set(yticks=[], ylabel="")
g.despine(bottom=True, left=True)

for ax in g.axes.flat:
    ax.axvline(x=df[value].median(), color='r', linestyle='-')


plt.show()


# In[31]:


def mad_outlier(adata, metric, nmads, upper_only = False):
    M = adata.obs[metric]
    
    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
    
    return (M > np.median(M) + nmads * mad(M))


# In[32]:


def pp(adata):
    adata = adata[adata.obs.pct_counts_mt < 25] #you can lower this based on the overal distribution of your dataset
    
    bool_vector = mad_outlier(adata, 'log1p_total_counts', 5) +\
            mad_outlier(adata, 'log1p_n_genes_by_counts', 5) +\
            mad_outlier(adata, 'pct_counts_in_top_20_genes', 5) +\
            mad_outlier(adata, 'pct_counts_mt', 3, upper_only = True)
    adata = adata[~bool_vector]

    adata.uns['cells_removed'] = sum(bool_vector)

    return adata


# In[33]:


adatas = [pp(ad) for ad in adatas]


# In[34]:


for adata in adatas:
    print(len(adata), adata.uns['cells_removed'])


# In[35]:


adatas


# In[37]:


df2 = pd.concat([x.obs for x in adatas])
df2 = df2.sort_values('Sample')


# In[38]:


#value = "pct_counts_mt"
value = "n_genes"
#value = 'pct_counts_in_top_20_genes'
#value = "log1p_total_counts"

sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

g = sns.FacetGrid(df2, row="Sample", hue="Sample", aspect=10, height=1, palette="tab20")

g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

g.map(plt.axhline, y=0, lw=2, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)

g.map(label, value)

g.figure.subplots_adjust(hspace=-0.2)

g.set_titles("")
g.set(yticks=[], ylabel="")
g.despine(bottom=True, left=True)

for ax in g.axes.flat:
    ax.axvline(x=df2[value].median(), color='r', linestyle='-')


plt.show()


# # Doublet removal

# In[39]:


import scanpy as sc
import scvi

# Iterar sobre las listas de AnnData en adatas
for i, adata in enumerate(adatas):
    print(f"Procesando el conjunto de datos {i+1}")
    
    # Paso 1: Preprocesamiento
    # Filtrar genes que no se expresan en al menos 10 células
    sc.pp.filter_genes(adata, min_cells=10)
    
    # Seleccionar los 2000 genes más variables
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False, flavor='seurat_v3')
    sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)
    # Verificación después del preprocesamiento
    print(f"Datos después del filtrado de genes y selección de genes variables: {adata.shape}")
    adata.X = adata.X.tocsr()
    # Paso 2: Configurar y entrenar el modelo scVI
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()

    # Paso 3: Aplicar SOLO para detectar dobletes
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()

    # Obtener las predicciones
    df = solo.predict()
    
    # Generar etiquetas binarias (predicción directa de doublet o singlet)
    df['prediction'] = solo.predict(soft=False)
    
    # Calcular la diferencia entre la probabilidad de 'doublet' y 'singlet'
    df['dif'] = df.doublet - df.singlet

    # Paso 4: Filtrar células dobletes con un umbral de dif > 0.9
    doublets = df[(df.prediction == 'doublet') & (df.dif > 0.9)]

    # Ver cuántas células están etiquetadas como dobletes
    print(f"Total de células etiquetadas como doblete con dif > 0.9: {doublets.shape[0]}")

    # Paso 5: Crear una nueva columna en adata.obs que indique si la célula es un doblete
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)

    # Paso 6: Filtrar las células que no son dobletes (i.e., mantener solo singletes)
    adata_filtrado = adata[~adata.obs['doublet']]

    # Actualizar la lista adatas con el objeto filtrado
    adatas[i] = adata_filtrado

    # Verificación después del filtrado de dobletes
    print(f"Datos restantes después de filtrar dobletes: {adata_filtrado.shape}")


# In[40]:


adata = sc.concat(adatas, join='outer')


# In[41]:


adata


# In[42]:


adata.write('/home/mdiaz/sc_liver_data/checkpoints/doublets_removed.h5ad')


# In[ ]:





# In[2]:


get_ipython().run_line_magic('cd', '/home/mdiaz/hackaton-csbig/ref_Data')


# In[3]:


get_ipython().system('wget https://datasets.cellxgene.cziscience.com/a48b7e8a-9db3-45f1-9729-c87738c0082f.h5ad')


# In[8]:


get_ipython().system('wget https://datasets.cellxgene.cziscience.com/6c7f190c-efc1-4d63-a272-26769dd1d1d1.h5ad')


# In[12]:


get_ipython().system('wget https://datasets.cellxgene.cziscience.com/8b1e0d38-c48a-4e60-84c4-744a1d8ad2a1.h5ad')


# In[15]:


get_ipython().system('wget https://datasets.cellxgene.cziscience.com/c7a60fef-4ea1-4415-87e3-65f2b5f3f89f.h5ad')


# In[16]:


get_ipython().system('ls')


# In[5]:


ref = sc.read_h5ad('a48b7e8a-9db3-45f1-9729-c87738c0082f.h5ad')


# In[6]:


ref


# In[7]:


ref.obs


# In[11]:


ref2 = sc.read_h5ad('6c7f190c-efc1-4d63-a272-26769dd1d1d1.h5ad')
ref2.obs


# In[14]:


ref3 = sc.read_h5ad('8b1e0d38-c48a-4e60-84c4-744a1d8ad2a1.h5ad')
ref3.obs


# In[18]:


ref4 = sc.read_h5ad('c7a60fef-4ea1-4415-87e3-65f2b5f3f89f.h5ad')
ref4


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




