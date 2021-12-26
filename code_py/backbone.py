import pandas as pd
import numpy as np
import sys 
sys.path.append('../')
from tqdm import tqdm
import networkx as nx
import markov_clustering as mc
import networkx as nx
import random
from tqdm import tqdm
from sklearn.model_selection import KFold
from code_py.DIAMOnD import *
from code_py.backbone import *
from sklearn.preprocessing import normalize
sys_path = "/Users/alessandroquattrociocchi/Documents/Courses /2/2.1/BI/final project/"
data_path = sys_path + "data/"

class Disease_Genes_Graph:

    def __init__(self,folder_path="/Users/alessandroquattrociocchi/Documents/Courses /2/2.1/BI/final project/"):
        self.folder_path = folder_path
        self.data_path   = folder_path + "data/"
        self.columns = ['Official Symbol Interactor A','Official Symbol Interactor B']
        super(Disease_Genes_Graph, self).__init__()
        
    
    def preprocessing_dataset(self,homo_sap=True,drop_duplicates=True,remove_self_loops=True):
        self.homo_sapiens_genes = pd.read_csv(self.data_path+'BIOGRID-ORGANISM-Homo_sapiens-4.4.204.tab3.txt', sep='\t', header=0,low_memory=False)
        if homo_sap:
            self.homo_sapiens_genes = self.homo_sapiens_genes[(self.homo_sapiens_genes["Organism Name Interactor A"]=='Homo sapiens') & (self.homo_sapiens_genes["Organism Name Interactor B"]=='Homo sapiens')]
            self.homo_sapiens_genes = self.homo_sapiens_genes[self.columns]
        if drop_duplicates:
            self.homo_sapiens_genes = self.homo_sapiens_genes.drop_duplicates()
        if remove_self_loops:
            self.homo_sapiens_genes = self.homo_sapiens_genes[(self.homo_sapiens_genes['Official Symbol Interactor A'] != self.homo_sapiens_genes['Official Symbol Interactor B'])]
        print("Number of putative genes:",self.homo_sapiens_genes.shape[0])
        return self.homo_sapiens_genes

    def query_disease_genes(self,diseaseId):
        #"C1510586"
        self.diseases = pd.read_csv(self.data_path+"curated_gene_disease_associations.tsv", sep='\t')
        self.disease_query =  self.diseases[self.diseases["diseaseId"]==diseaseId]
        self.disease_list = list(self.disease_query['geneSymbol'])
        print("Found " + str(len(self.disease_list)) + ' disease genes in ' + str(self.disease_query['diseaseName'].values[0]))
        return self.disease_query,self.disease_list

    def LCC_to_adj(self,dataframe):
        self.putative_genes_graph = nx.from_pandas_edgelist(dataframe, source = "Official Symbol Interactor A", target = "Official Symbol Interactor B", 
                              create_using=nx.Graph())

        #finding the connected components
        self.conn_comp = list(nx.connected_components(self.putative_genes_graph))
        #len of the connected component 
        self.conn_comp_len = [len(c) for c in sorted(self.conn_comp, key=len, reverse=True)]
        print("# of connected components:", len(self.conn_comp))
        #finding the largest connected component 
        self.LCC = max(self.conn_comp, key=len)
        print(len(self.LCC)) #LCC is a dict of nodes 
        #creating a subgraph with the largest connected component
        self.LCC_sub_graph = self.putative_genes_graph.subgraph(self.LCC).copy()               
        print(nx.info(self.LCC_sub_graph))
        #converting subgraph into the adj matrix 
        self.LCC_sub_graph_adj = nx.adjacency_matrix(self.LCC_sub_graph)
        return self.LCC_sub_graph,self.LCC_sub_graph_adj , self.putative_genes_graph.number_of_nodes(), self.putative_genes_graph.number_of_edges()

    @staticmethod
    def MCL(adj_mat,inflation):
        cluster_data = {}
        count = 0
        result = mc.run_mcl(adj_mat, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        count += 1
        return Q
        
    @staticmethod
    def KFold_CV(list,n_folds=10,shuffle_flag=True):
        kf = KFold(n_splits=n_folds, shuffle=shuffle_flag)
        X = np.array(list)
        X_dataset_cv = []
        Y_dataset_cv = []
        for train_index, val_index in kf.split(X):
            X_dataset_cv.append(X[train_index].tolist())
            Y_dataset_cv.append(X[val_index].tolist())
        return X_dataset_cv,Y_dataset_cv
    
    @staticmethod
    def RWR(LCC_sub_graph,restart_prob = 0.5,conv_treshold=0.000001,max_print_items= 10):
        source = list(LCC_sub_graph.nodes())
        p_0 = [0]*LCC_sub_graph.number_of_nodes()
        for source_id in source:
            source_index = list(LCC_sub_graph.nodes()).index(source_id)
            p_0[source_index] = 1 / float(len(source))
        p_0 = np.array(p_0)
        og_matrix = nx.to_numpy_matrix(LCC_sub_graph)
        og_matrix = normalize(og_matrix, norm='l1', axis=0)
        
        def calculate_next_p(og_matrix, p_t, p_0):
            """ Calculate the next probability vector. """
            epsilon = np.squeeze(np.asarray(np.dot(og_matrix, p_t)))
            no_restart = epsilon * (1 -restart_prob)
            restart = p_0 * restart_prob
            return np.add(no_restart, restart)
        
        # set up the starting probability vector
        diff_norm = 100
        # this needs to be a deep copy, since we're reusing p_0 later
        p_t = np.copy(p_0)

        while (diff_norm > conv_treshold):
            # first, calculate p^(t + 1) from p^(t)
            p_t_1 = calculate_next_p(og_matrix,p_t, p_0)

            # calculate L1 norm of difference between p^(t + 1) and p^(t),
            # for checking the convergence condition
            diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)

            # then, set p^(t) = p^(t + 1), and loop again if necessary
            # no deep copy necessary here, we're just renaming p
            p_t = p_t_1
        gene_probs =zip(list(LCC_sub_graph.nodes()), p_t.tolist())
        s = sorted(gene_probs, key=lambda x: x[1], reverse=True)
        for i in range(0,max_print_items):
            print(s[i][0],s[i][1])