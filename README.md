# paper deep learning vs gblup
Datasets and scripts used for the article "Deep learning vs GBLUP for whole-genome predictions: relative efficiency with additive and non-additive continuous phenotypes"

- `deep_learning_model_run.ipynb` is the `Python` notebook with the code to run the deep learning model developed in the paper. This is `tensorflow`/`keras` implementation, which can be run on [Google Colaboratory](https://colab.research.google.com/) (or a similar Jupyter Notebook engine)
- the folder `support_scripts` contains all `Python` dependencies needed to run the DL model (download and prepare the data, build the model, parse and collect results)
- the data are stored in this repository (simulated phenotypes) and on [zenodo.org](https://zenodo.org/record/6602439#.YpofCHVBxhE) (kinship matrices)
- the `R` script to simulate the phenotypes with different relative contributions of additive, dominant and epistatic genetic effects is also included (in `support_scripts`)
