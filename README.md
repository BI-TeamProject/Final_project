# Comparative assessment of disease gene prediction algorithms

<p align="center">
<img src= https://github.com/BI-TeamProject/Final_project/blob/main/Further%20Material/git_images/sapienza_logo.jpg
 width="100"/>
 </p>
  
  <p align="center">
  <b>BIOINFORMATICS AND NETWORK MEDICINE a.y 2021-2022<br />
La Sapienza University of Rome <br />
MSC IN DATA SCIENCE<b>  <br />
Alessandro Quattrociocchi, Tansel Simsek<b> <br />
</p>
  
  

## Abstract
In this report, we will interview different methodologies to investigate protein-protein interactions and its' role in the disease gene association. The aim is to find the most useful algorithm to extract correct genes to trace disease existence. We consider 5 different diseases, 5 different algorithms which are MCL (Markov Cluster), DIAMOnD (a disease module detection), DiaBLE, heat diffusion with Cytoscape bioinformatics software and Random Walk with Restart, and some of evaluation metrics such as precision, recall, F1 score and nDCG(normalized discounted cumulative gain). Also, we tried 6 additional diseases which all the results can be found in the additional pdf. As a final result, we investigated the most suitable disease as Autism Spectrum Disorder and the most powerful algorithm as DiaBLE with further enrichment analysis.

<p align="center">
<img src=https://github.com/BI-TeamProject/Final_project/blob/main/Further%20Material/git_images/ppi_image.jpg width="500"/ >
</p>
   
### Folder Tree
```bash
├── code_py
│   ├── DIAMOnD.py
│   ├── __init__.py
│   └── backbone.py
├── main.ipynb
└── outputs
    ├── MLC_modularity
    ├── pkl_datasets
    │   ├── C0020796.pkl
    │   ├── C0079744.pkl
    │   ├── C0086565.pkl
    │   ├── C0205644.pkl
    │   ├── C0238198.pkl
    │   ├── C0860207.pkl
    │   ├── C1510586.pkl
    │   ├── C1959583.pkl
    │   ├── C3714756.pkl
    │   ├── C4316881.pkl
    │   └── C4505456.pkl
    └── results_table
        ├── C0020796.html
        ├── C0079744.html
        ├── C0086565.html
        ├── C0205644.html
        ├── C0238198.html
        ├── C0860207.html
        ├── C1510586.html
        ├── C1959583.html
        ├── C3714756.html
        ├── C4316881.html
        └── C4505456.html

4 directories, 27 files
```

   #### Image Credits 
   https://deepmind.com/blog/article/putting-the-power-of-alphafold-into-the-worlds-hands
