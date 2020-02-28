import warnings
import sys
import os
import torch
import matplotlib
warnings.filterwarnings('ignore')

sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.inference import UnsupervisedTrainer, AlternateSemiSupervisedTrainer, SemiSupervisedTrainer
from scvi.models import VAE, SCANVI
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def get_scvi_posterior(data, model_file, retrain=False, seed=0, n_epochs=150,
                       use_batches=True, use_cuda=True, lr=1e-3):
    vae = VAE(data.nb_genes, n_batch=data.n_batches * use_batches,
              n_layers=2, n_latent=30, dispersion='gene')
    trainer = UnsupervisedTrainer(
        vae,
        data,
        train_size=0.99,
        use_cuda=use_cuda,
        frequency=5,
        data_loader_kwargs={"pin_memory": False},
        n_iter_kl_warmup = 1600
    )

    torch.manual_seed(seed)
    if retrain == True:
        trainer.train(n_epochs=n_epochs, lr=lr)
        torch.save(trainer.model.state_dict(), model_file)
    else:
        if os.path.isfile(model_file):
            trainer.model.load_state_dict(torch.load(model_file))
            trainer.model.eval()
        else:
            trainer.train(n_epochs=n_epochs, lr=lr)
            torch.save(trainer.model.state_dict(), model_file)

    elbo_train_set = trainer.history["elbo_train_set"]
    elbo_test_set = trainer.history["elbo_test_set"]
    x = np.linspace(0, 150, (len(elbo_train_set)))
    plt.plot(x[2:], elbo_train_set[2:])
    plt.plot(x[2:], elbo_test_set[2:])
    plt.savefig(model_file + '.ll.pdf')
    posterior = trainer.create_posterior(
        trainer.model, data, indices=np.arange(len(data))
    ).sequential()

    return posterior


def balance_n_labelled(l, labelled, nlabels=30):
    balanced_labelled = []
    for i in np.unique(l[labelled]):
        idx = np.where(l == i)[0]
        if len(idx) >= nlabels:
            subset = np.random.choice(idx, nlabels, replace=False)
        else:
            subset = idx
        balanced_labelled.append(subset)
    balanced_labelled = np.concatenate(balanced_labelled)
    return (balanced_labelled)


def scanvi_pred(data, modelfile, scanvi_modelfile, balance=True, alternate=False, nlabels=30,
                retrain=False, n_epochs=15, forward_only=False):
    if forward_only:
        print('Only predict based on pretrained model')
    else:
        if 'unassigned' in data.cell_types:
            unlabelled_idx = list(data.cell_types).index('unassigned')
            labelled = np.where(data.labels.ravel()!=unlabelled_idx)[0]
        else:
            labelled = np.arange(len(data.labels))

        # balance number of labelled cells from each cell type
        if balance is True:
            labelled = balance_n_labelled(data.labels.ravel(), labelled, nlabels=nlabels)
        labelled = np.random.choice(labelled, len(labelled), replace=False)

        unlabelled = [x for x in np.arange(len(data.labels.ravel())) if x not in labelled]
        unlabelled = np.random.choice(unlabelled, len(unlabelled), replace=False)

    scanvi = SCANVI(data.nb_genes, data.n_batches, data.n_labels, n_layers=2, n_latent=30,
                    symmetric_kl=True)

    if alternate is False:
        trainer_scanvi = SemiSupervisedTrainer(scanvi, data,
                                               n_epochs_classifier=100,
                                               lr_classification=5 * 1e-3, seed=1,
                                               n_epochs_kl_warmup=1)
    else:
        trainer_scanvi = AlternateSemiSupervisedTrainer(scanvi, data,
                                                        n_epochs_classifier=100,
                                                        lr_classification=5 * 1e-3,
                                                        n_epochs_kl_0warmup=1)

    if forward_only:
        trainer_scanvi.model.load_state_dict(torch.load(scanvi_modelfile))
        trainer_scanvi.model.eval()
    else:
        trainer_scanvi.model.load_state_dict(torch.load(modelfile), strict=False)
        trainer_scanvi.model.eval()
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=labelled)
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=unlabelled)
        if retrain:
            trainer_scanvi.train(n_epochs=n_epochs)
            torch.save(trainer_scanvi.model.state_dict(), scanvi_modelfile)
        else:
            if os.path.isfile(scanvi_modelfile):
                trainer_scanvi.model.load_state_dict(torch.load(scanvi_modelfile))
                trainer_scanvi.model.eval()
            else:
                trainer_scanvi.train(n_epochs=n_epochs)
                torch.save(trainer_scanvi.model.state_dict(), scanvi_modelfile)

    full = trainer_scanvi.create_posterior(trainer_scanvi.model, data, indices=np.arange(len(data)))
    _, pred = full.sequential().compute_predictions()
    pred_celltype = [data.cell_types[i] for i in pred]

    return full, pred_celltype


def DEbyCompartment(raw, full, label, split, split_subset, filename):
    pred_label = np.asarray(raw.obs[label].values)
    pred_celltype, pred_label = np.unique(pred_label, return_inverse=True)
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for compartment in split_subset:
        de_res, de_cluster = full.one_vs_all_degenes(cell_labels=pred_label,
                                                     subset=raw.obs[split].values == compartment)

        for i, ct in enumerate(pred_celltype[de_cluster]):
            x = de_res[i]
            filt = np.logical_and(x['bayes_factor'].values > 1.3, x['raw_mean1'].values > 0.1)
            res = x.loc[filt]
            res.to_excel(writer, sheet_name=ct)
            res['celltype'] = ct

    writer.save()


def check_missing_celltypes(adata, pred):
    missing = [x for x in np.unique(adata.obs['input_ann2'])
               if x not in np.unique(adata.obs[pred]) and x != 'unassigned']
    missing10x = [x for x in np.unique(adata.obs['input_ann2'])
                  if x not in np.unique(np.asarray(adata.obs[pred])[adata.obs['batch'] == '0']) and x != 'unassigned']
    missingss2 = [x for x in np.unique(adata.obs['input_ann2'])
                  if x not in np.unique(np.asarray(adata.obs[pred])[adata.obs['batch'] == '1']) and x != 'unassigned']
    res = "%s missing,\n %s missing from 10x,\n %s missing from Smartseq2" % (
        ", ".join(missing), ", ".join(missing10x), ", ".join(missingss2))
    return res


def clustercomp(List):
    List = list(List)
    x, k = np.unique(List,return_counts=True)
    if np.max(k)/np.sum(k)>0.9:
        return max(set(List), key = List.count)
    else:
        return('mixed')


def HierarchicalFrequency(celltypes, freq, co):
    ont = {}
    for x in co.nodes.keys():
        for celltype in celltypes:
            if co.nodes[x]['name'] is celltype:
                ont[celltype] = x
    old_freq = dict(zip(celltypes, freq))
    new_freq = dict(zip(celltypes, freq))
    for x in celltypes:
        children = ancestors(co, ont[x])
        for y in celltypes:
            if ont[y] in children:
                new_freq[x] += old_freq[y]
    new_freq = np.asarray([new_freq[x] for x in celltypes])
    return new_freq


from sklearn.neighbors import NearestNeighbors
import scipy


def entropy(hist_data):
    n_batches = len(np.unique(hist_data))
    if n_batches > 2:
        raise ValueError("Should be only two clusters for this metric")
    frequency = np.mean(hist_data == 1)
    if frequency == 0 or frequency == 1:
        return 0
    return -frequency * np.log(frequency) - (1 - frequency) * np.log(1 - frequency)


def entropy_batch_mixing(latent, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100):
    n_neighbors = min(n_neighbors, len(latent) - 1)
    nne = NearestNeighbors(n_neighbors=1 + n_neighbors, n_jobs=8)
    nne.fit(latent)
    kmatrix = nne.kneighbors_graph(latent) - scipy.sparse.identity(
        latent.shape[0]
    )
    score = 0
    for t in range(n_pools):
        indices = np.random.choice(
            np.arange(kmatrix.shape[0]), size=n_samples_per_pool
        )
        pool = 0
        for i in range(n_samples_per_pool):
            nonzeros = kmatrix[indices].nonzero()
            sampled_cell = nonzeros[0]
            all_neighbors = nonzeros[1]
            neighbors = all_neighbors[sampled_cell==i]
            pool += entropy(batches[neighbors])
        score += pool/n_samples_per_pool

    return score / float(n_pools)


def entropy_by_category(adata, latent_name, label_col, organ_col, n_neighbors=100):
    organs = adata.obs[organ_col].values
    celltypes = adata.obs[label_col].values
    entropy = pd.DataFrame(index = np.unique(celltypes))
    latent = adata.obsm[latent_name]
    for o in np.unique(organs):
        entropy[o+'_organ_mixing'] = 0
        for c in np.unique(celltypes[organs==o]):
            res = entropy_batch_mixing(latent[celltypes==c,:],(organs[celltypes==c]==o).astype(int), n_neighbors=n_neighbors)
            print(o,c, res)
            entropy.loc[c, o+'_organ_mixing'] = res
    return(entropy)
