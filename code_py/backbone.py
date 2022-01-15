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
from sklearn.preprocessing import normalize
from operator import itemgetter
from scipy.stats import hypergeom
import pickle
import statistics 
import math
import statistics 


class Human_Genes_Graph_Analysis:

    def __init__(self,folder_path,disease_ID):
        self.folder_path = folder_path
        self.data_path   = folder_path + "data/"
        self.disease_ID = disease_ID
        super(Human_Genes_Graph_Analysis, self).__init__()
    
    # ============================  PREPROCESSING 
    def preprocessing_dataset(self, homo_sap=True, drop_duplicates=True, remove_self_loops=True, write_txt = True):
        """
        filtering dataset
        :params homo_sap: bool -> filtering the dataset accordingly homo sapiens genes
        :params drop_duplicates: bool -> removes the dusplicates in the dataset
        :params remove_self_loops: bool -> removes the self loops from the ppi
        :write_txt: bool -> writes output txt file
        """
        self.homo_sapiens_genes = pd.read_csv(self.data_path+'BIOGRID-ORGANISM-Homo_sapiens-4.4.204.tab3.txt', sep='\t', header=0,low_memory=False)
        if homo_sap:
            self.homo_sapiens_genes = self.homo_sapiens_genes[(self.homo_sapiens_genes["Experimental System Type"]=='physical')]
            self.homo_sapiens_genes = self.homo_sapiens_genes[(self.homo_sapiens_genes["Organism ID Interactor A"]==9606) & (self.homo_sapiens_genes["Organism ID Interactor B"]==9606)]
            self.trial = self.homo_sapiens_genes
        if write_txt:
            self.homo_sapiens_genes[['Official Symbol Interactor A', 'Official Symbol Interactor B']].to_csv(self.folder_path +'data/Biogrid_4.4.204.txt', header=None, index=None, sep=' ', mode='a')
        if drop_duplicates:
            self.homo_sapiens_genes = self.homo_sapiens_genes.drop_duplicates()
        if remove_self_loops:
            self.homo_sapiens_genes = self.homo_sapiens_genes[(self.homo_sapiens_genes['Official Symbol Interactor A'] != self.homo_sapiens_genes['Official Symbol Interactor B'])]
        return self.homo_sapiens_genes
    
    # ============================  QUERY DISEASE GENES 
    def query_disease_genes(self):
        """
        Filter on the diseases curated dataset according to the input disease query
        """
        self.diseases = pd.read_csv(self.data_path+"curated_gene_disease_associations.tsv", sep='\t')
        self.disease_query =  self.diseases[self.diseases["diseaseId"]==self.disease_ID]
        self.disease_list = list(self.disease_query['geneSymbol'])
        print("Found " + str(len(self.disease_list)) + ' disease genes in ' + str(self.disease_query['diseaseName'].values[0]))
        return self.disease_query,self.disease_list

    # ============================ QUERY DISEASE GENES EXTENDED
    def query_disease_genes_extendend(self):
        """
        Filter on the diseases all genes associatiations dataset according to the input disease query
        """
        self.diseases_ex = pd.read_csv(self.data_path+"all_gene_disease_associations.tsv", sep='\t')
        self.diseases_ex_filtered = self.diseases_ex[self.diseases_ex['diseaseId']==self.disease_ID]
        return self.diseases_ex_filtered,list(self.diseases_ex_filtered['geneSymbol'].values)

    # ============================  LARGE CONNECTED COMPONENTS TO ADJ 
    def LCC_to_adj(self,dataframe):
        """
        Converting tha dataframe as graph with targets and sources, returns the LCC sub graph, adjacency matrix, # nodes and #edges
        :param: datasframe
        """
        self.putative_genes_graph = nx.from_pandas_edgelist(dataframe, source = "Official Symbol Interactor A", target = "Official Symbol Interactor B", 
                              create_using=nx.Graph())
        #finding the connected components
        self.conn_comp = list(nx.connected_components(self.putative_genes_graph))
        #len of the connected component 
        self.conn_comp_len = [len(c) for c in sorted(self.conn_comp, key=len, reverse=True)]
        #finding the largest connected component 
        self.LCC = max(self.conn_comp, key=len)
        print(len(self.LCC)) #LCC is a dict of nodes 
        #creating a subgraph with the largest connected component
        self.LCC_sub_graph = self.putative_genes_graph.subgraph(self.LCC).copy()               
        print(nx.info(self.LCC_sub_graph))
        #converting subgraph into the adj matrix 
        self.LCC_sub_graph_adj = nx.adjacency_matrix(self.LCC_sub_graph)
        return self.LCC_sub_graph,self.LCC_sub_graph_adj , self.putative_genes_graph.number_of_nodes(), self.putative_genes_graph.number_of_edges()

    # ============================  K-FOLD CROSS VALIDATION
    def KFold_CV(self, list, n_folds=5, shuffle_flag=True):
        """
        K-Fold cross validation
        :param list: -> input list 
        :param n_folds: int -> default=5
        :param shuffle_flag: bool -> default= True
        :return X,Y: lists with folds
        """
        kf = KFold(n_splits=n_folds, shuffle=shuffle_flag,random_state=1234567)
        X = np.array(list)
        X_dataset_cv = []
        Y_dataset_cv = []
        for train_index, val_index in kf.split(X):
            X_dataset_cv.append(X[train_index].tolist())
            Y_dataset_cv.append(X[val_index].tolist())
        return X_dataset_cv,Y_dataset_cv

    # ============================  MCL ALGORITHM 
    @staticmethod
    def MCL(adj_mat,inflation):
        """
        MCL modularity computation
        :params adj_mat: array -> adjacency matrix of the graph
        :param inflation: int
        :return Q:float -> modularity value 
        """
        count = 0
        result = mc.run_mcl(adj_mat, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        count += 1
        return Q
    # ============================ MCL CLUSTERS EVALUATION
    @staticmethod
    def MLC_eval(disease_graph,disease_genes_list,clusters,tol=0):
        """
        MCL clusters evaluation 
        :param disease_graph: graph object 
        :param disease_genes_list: input disease list
        :param clusters: list
        :return pmf_cluster_dict_avg: dict with index cluster as key, p-values as values
        :return enriched_genes_list: list -> enriched clusters list
        :return enriched_cluster_index: list -> list of the index of the enriched clusters
        """
        pmf_cluster_dict = {}
        for folds_index in range(0,5):
            print("================================================")
            print("Fold number: ",folds_index )
            for cluster_index in range(0,len(clusters)): 
                ds_genese_as_string = set(itemgetter(*clusters[cluster_index])(list(disease_graph.nodes)))
                intersect_cluster_size = len(set(disease_genes_list[folds_index]).intersection(ds_genese_as_string))

                if intersect_cluster_size > 2:
                    
                    [N, K, n] = [19618, len(disease_genes_list[folds_index]), len(ds_genese_as_string)]
                    rv = hypergeom(N, K, len(clusters[cluster_index]))
                    pmf_cluster = rv.pmf(intersect_cluster_size)
                    print(str(intersect_cluster_size) + " disease genes " + "in" + " cluster " + str(cluster_index) + " --> " + str(round(pmf_cluster,6)))
                    if cluster_index not in list(pmf_cluster_dict.keys()):
                        pmf_cluster_dict[cluster_index] = [pmf_cluster]
                    elif cluster_index in list(pmf_cluster_dict.keys()):
                        pmf_cluster_dict[cluster_index].append(pmf_cluster)

        pmf_cluster_dict_avg = {}
        for k,v in pmf_cluster_dict.items():
            # v is the list of grades for student k
            pmf_cluster_dict_avg[k] = sum(v)/ float(len(v))

        enriched_cluster_index = [k for k,v in pmf_cluster_dict_avg.items() if float(v) <= 0.05+tol]
        print("================================================")
        print("The index of the enriched cluster found using MLC is: ",enriched_cluster_index)
        enriched_genes_list = [] 
        for c in enriched_cluster_index:

            enriched_genes_list.append(set(itemgetter(*clusters[c])(list(disease_graph.nodes))))
        return pmf_cluster_dict_avg, enriched_genes_list, enriched_cluster_index
    
    # ============================ MCL METRICS EVALUATION
    @staticmethod
    def MCL_evaluation_metrics(disease_graph,dg_list,hs_disease_genes,clusters,enriched_clusters_list):
        """
        MCL evaluation 
        :param disease_graph: graph object
        :param dg_list: list-> diseases genes
        :param hs_disease_genes: list
        :param clusters: list  -> list of clusters
        :param enriched_clusters_list: -> enriched clsuter list
        """
        genes_in_all_clusters = []
        TP = 0 
        FP = 0 
        intersect_clusters = []
        for fold in dg_list:
            for e in enriched_clusters_list:
                ds_genese_as_string = set(itemgetter(*clusters[e])(list(disease_graph.nodes)))
                genes_in_all_clusters += list(ds_genese_as_string)
                intersect_clusters.extend(list(set(fold).intersection(ds_genese_as_string)))
            
        #number of probe set genes in all enriched clusters
        TP = len(set(intersect_clusters))
        #number of genes in all enriched cluster which are not seed genes
        FP = len(set(clusters))-TP
        #number of probe seed genes not present in any enriched clusters
        FN = len(set(disease_graph.nodes))-len(set(genes_in_all_clusters).intersection(set(hs_disease_genes)))
        print("TP: " + str(TP) + " --- " + "FP: " +str(FP) + " --- " + "FN: " +str(FN))
        precision=(TP/(TP+FP))
        recall   = ((TP)/(TP+FN))
        f1_score = ((2*precision*recall)/(precision+recall))
        print("Precision: " + str(round(precision*100,2)) + " --- " + "Recall: " +str(round(recall*100,2)) + " --- " + "F1 Score: " +str(round(f1_score*100,2)))

    # ============================ RANDOMD WALK WITH RESTART 
    def RWR(self, LCC_sub_graph,disease_genes,restart_prob = 0.75,conv_treshold=1e-6,max_print_items= 10):
        """
        Random Walk With Restart Algorithm
        :param LCC_sub_graph: graph object 
        :param disease_genes: list 
        :param restart_prob:int -> restarting prob value
        :param conv_treshold: float32 -> convergence threshold 
        :param max_print_items: int -> max printed element in output
        :return s: list -> ranking of putative genes
        """
        source = set(disease_genes).intersection(set(LCC_sub_graph.nodes()))
        p_0 = [0]*len(list(LCC_sub_graph.nodes()))
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
        genens_ranked = set([list(i)[0] for i in s]) - set(disease_genes)

        return s,list(genens_ranked)

    # ============================ METRICS EVALUATION  ENGINE

    def ndcg(self, predicted_list, true_list, n):
        """
        normalized discounted cumulative gain calculation at top n position
        :param predicted_list: list -> contains the predicted disease genes/putative genes 
        :param true_list: list -> ground truth of the disease that contain true disease genes
        :param n:int -> top n position 
        :return dcg/idcg: number -> nDCG score
        """
        dcg = 0
        idcg = 0
        for i in range(n):
            if predicted_list[i] in true_list:
                dcg += 1/math.log2(i+ 2)
            idcg += 1/math.log2(i+ 2)
        return dcg/idcg
    
    def return_pre_recall(self, predicted_nodes, ds_genes_train, ds_genes_test, n):
        """
        precision and recall calculation at top n position
        :param predicted_nodes: list -> contains the predicted disease genes/putative genes 
        :param ds_genes_train: list -> train dataset of the seed genes
        :param ds_genes_test: list -> test dataset of the seed genes
        :param n:int -> top n position 
        :return precision and recall values for top 50, n/10, n/4, n/2 and n positions
        """
        TP_50   = len(set(predicted_nodes[:50]).intersection(set(ds_genes_test)))
        TP_n_10 = len(set(predicted_nodes[:n//10]).intersection(set(ds_genes_test)))
        TP_n_4  = len(set(predicted_nodes[:n//4]).intersection(set(ds_genes_test)))
        TP_n_2  = len(set(predicted_nodes[:n//2]).intersection(set(ds_genes_test)))
        TP_n = len(set(predicted_nodes[:min(n, len(predicted_nodes))]).intersection(set(ds_genes_test)))

        FP_50   = (len(set(predicted_nodes[:50]))-TP_50)
        FP_n_10 = (len(set(predicted_nodes[:n//10])) - TP_n_10)
        FP_n_4  =  (len(set(predicted_nodes[:n//4]))  - TP_n_4)
        FP_n_2  =  (len(set(predicted_nodes[:n//2]))  - TP_n_2)
        FP_n    =  (len(set(predicted_nodes[:min(n, len(predicted_nodes))])) - TP_n)
        
        FN_50 = len(ds_genes_test) - TP_50
        FN_n_10 = len(ds_genes_test) - TP_n_10
        FN_n_4 = len(ds_genes_test) - TP_n_4
        FN_n_2 = len(ds_genes_test) - TP_n_2
        FN_n = len(ds_genes_test) - TP_n
        
        return TP_50/(TP_50+FP_50), TP_n_10/(TP_n_10+FP_n_10), TP_n_4/(TP_n_4+FP_n_4), TP_n_2/(TP_n_2+FP_n_2), TP_n/(TP_n+FP_n), TP_50/(TP_50+FN_50), TP_n_10/(TP_n_10+FN_n_10), TP_n_4/(TP_n_4+FN_n_4), TP_n_2/(TP_n_2+FN_n_2), TP_n/(TP_n+FN_n) 

    def return_metrics(self, method, pgenes_sub_graph, hs_disease_genes, ds_genes_train, ds_genes_test, print_flag=True,save_dataframe=True, extended_val = False):
        """
        all metrics calculation for the ranked prediction
        :param method: string -> ["DIAMOnD", "DiaBLE", "Cytoscape", "RWR"] 
        :param hs_disease_genes: list -> contains the true disease genes
        :param ds_genes_train: list -> train dataset of the seed genes
        :param ds_genes_test: list -> test dataset of the seed genes
        :param extended_val:boolean -> True if extended validation part will be, False otherwise
        :return: dataframe -> contains precision, recall, F1 score and nDCG values for top 50, n/10, n/4, n/2 and n positions for a given disease
        """    
        precision_50,precision_n_10,precision_n_4,precision_n_2,precision_n = [],[],[],[],[]
        recall_50, recall_n_10, recall_n_4, recall_n_2, recall_n = [],[],[],[],[]
        f1_score_50, f1_score_n_10, f1_score_n_4, f1_score_n_2, f1_score_n = [],[],[],[],[]
        ndcg_50, ndcg_n_10, ndcg_n_4, ndcg_n_2, ndcg_n = [],[],[],[],[]
        n = len(hs_disease_genes)

        if extended_val:
            _, ev_disease_genes = self.query_disease_genes_extendend()
            n_ev = len(ev_disease_genes)
            ev_genes_train, ev_genes_test = self.KFold_CV(ev_disease_genes)
            
        predicted_nodes_all = []
        if method == "DIAMOnD":
            for i in range(5):
                if extended_val:
                    added_nodes, predicted_nodes = DIAMOnD(G_original=pgenes_sub_graph,
                                seed_genes=ds_genes_train[i],
                                max_number_of_added_nodes=n_ev,alpha=1,DiaBLE=False)
                    TP_list = set(predicted_nodes).intersection(set(ds_genes_test[i]))
                    FP_predicted_nodes = [i for i in predicted_nodes if i not in TP_list]
                    predicted_nodes_all.append(FP_predicted_nodes)
                else:
                    added_nodes, predicted_nodes = DIAMOnD(G_original=pgenes_sub_graph,
                                seed_genes=ds_genes_train[i],
                                max_number_of_added_nodes=n,alpha=1,DiaBLE=False)
                    predicted_nodes_all.append(predicted_nodes)
        elif method == "DiaBLE":
            for i in range(5):
                if extended_val:
                    added_nodes, predicted_nodes = DIAMOnD(G_original=pgenes_sub_graph,
                        seed_genes=ds_genes_train[i],
                        max_number_of_added_nodes=n_ev,alpha=1,
                        DiaBLE=True)
                    TP_list = set(predicted_nodes).intersection(set(ds_genes_test[i]))
                    FP_predicted_nodes = [i for i in predicted_nodes if i not in TP_list]
                    predicted_nodes_all.append(FP_predicted_nodes)
                else:
                    added_nodes, predicted_nodes = DIAMOnD(G_original=pgenes_sub_graph,
                        seed_genes=ds_genes_train[i],
                        max_number_of_added_nodes=n,alpha=1,
                        DiaBLE=True)
                    predicted_nodes_all.append(predicted_nodes)
        elif method == "Cytoscape":
            for i in range(5):
                df = pd.read_csv('cytoscape/'+str(self.disease_ID)+'/diff_'+str(i)+'.csv')
                df = df[df.selected == False]
                df.sort_values(by=['diffusion_output_rank'], inplace=True)
                predicted_nodes = list(df['name'])
                if extended_val:
                    TP_list = set(predicted_nodes).intersection(set(ds_genes_test[i]))
                    FP_predicted_nodes = [i for i in predicted_nodes if i not in TP_list]
                    predicted_nodes_all.append(FP_predicted_nodes)
                else:
                    predicted_nodes_all.append(predicted_nodes)
        elif method == "RWR":
            for i in range(5):
                rwr_enriched_genes,_ = self.RWR(pgenes_sub_graph,ds_genes_train[i],max_print_items=0)
                predicted_nodes = [list(i)[0] for i in rwr_enriched_genes]
                if extended_val:
                    TP_list = set(predicted_nodes).intersection(set(ds_genes_test[i]))
                    FP_predicted_nodes = [i for i in predicted_nodes if i not in TP_list]
                    predicted_nodes_all.append(FP_predicted_nodes)
                else:
                    predicted_nodes_all.append(predicted_nodes)
        for i in range(len(predicted_nodes_all)):
            predicted_nodes = predicted_nodes_all[i]
            if extended_val:
                pre_50, pre_n_10, pre_n_4, pre_n_2, pre_n, rec_50, rec_n_10, rec_n_4, rec_n_2, rec_n = self.return_pre_recall(predicted_nodes, ev_genes_train[i], ev_genes_test[i], n_ev)
            else:
                pre_50, pre_n_10, pre_n_4, pre_n_2, pre_n, rec_50, rec_n_10, rec_n_4, rec_n_2, rec_n = self.return_pre_recall(predicted_nodes, ds_genes_train[i], ds_genes_test[i], n)
            precision_50.append(pre_50); precision_n_10.append(pre_n_10); precision_n_4.append(pre_n_4); precision_n_2.append(pre_n_2); precision_n.append(pre_n)
            recall_50.append(rec_50); recall_n_10.append(rec_n_10); recall_n_4.append(rec_n_4); recall_n_2.append(rec_n_2); recall_n.append(rec_n)
            if extended_val:
                ndcg_50.append(self.ndcg(predicted_nodes, ev_genes_test[i], 50))
                ndcg_n_10.append(self.ndcg(predicted_nodes, ev_genes_test[i], n_ev//10))
                ndcg_n_4.append(self.ndcg(predicted_nodes, ev_genes_test[i], n_ev//4))
                ndcg_n_2.append(self.ndcg(predicted_nodes, ev_genes_test[i], n_ev//2))
                ndcg_n.append(self.ndcg(predicted_nodes, ev_genes_test[i], min(n_ev, len(predicted_nodes))))
            else:
                ndcg_50.append(self.ndcg(predicted_nodes, ds_genes_test[i], 50))
                ndcg_n_10.append(self.ndcg(predicted_nodes, ds_genes_test[i], n//10))
                ndcg_n_4.append(self.ndcg(predicted_nodes, ds_genes_test[i], n//4))
                ndcg_n_2.append(self.ndcg(predicted_nodes, ds_genes_test[i], n//2))
                ndcg_n.append(self.ndcg(predicted_nodes, ds_genes_test[i], min(n, len(predicted_nodes))))
            try:
                f1_score_50.append((2*pre_50*rec_50)/(pre_50+rec_50))
            except:
                pass
            try:
                f1_score_n_10.append((2*pre_n_10*rec_n_10)/(pre_n_10+rec_n_10))
            except:
                pass
            try:
                f1_score_n_4.append((2*pre_n_4*rec_n_4)/(pre_n_4+rec_n_4))
            except:
                pass
            try:
                f1_score_n_2.append((2*pre_n_2*rec_n_2)/(pre_n_2+rec_n_2))
            except:
                pass
            try:
                f1_score_n.append((2*pre_n*rec_n)/(pre_n+rec_n))
            except:
                pass 
            
        if save_dataframe: 
  
            df = self.read_pkl_dataset(self.disease_ID)
            col_to_add  = method
            if extended_val: 
                col_to_add = col_to_add + " Ext"
                print(col_to_add)
            
            try:
                f1_str_50 = str(round(statistics.mean(f1_score_50)*100,2)) + " ± " +str(round(statistics.stdev(f1_score_50)*100,2))
            except:
                f1_str_50 = "0.0 ± 0.0"

            try:
                f1_str_n_10 = str(round(statistics.mean(f1_score_n_10)*100,2)) + " ± " +str(round(statistics.stdev(f1_score_n_10)*100,2))
            except:
                f1_str_n_10 = "0.0 ± 0.0"

            try:
                f1_str_n_4 = str(round(statistics.mean(f1_score_n_4)*100,2)) + " ± " +str(round(statistics.stdev(f1_score_n_4)*100,2))
            except:
                f1_str_n_4 = "0.0 ± 0.0" 

            try:
                f1_str_n_2 = str(round(statistics.mean(f1_score_n_2)*100,2)) + " ± " +str(round(statistics.stdev(f1_score_n_2)*100,2))
            except:
                f1_str_n_2 = "0.0 ± 0.0" 
            try: 
                f1_str_n = str(round(statistics.mean(f1_score_n)*100,2)) + " ± " +str(round(statistics.stdev(f1_score_n)*100,2))
            except:
                f1_str_n = "0.0 ± 0.0" 

            df = df.assign(**{col_to_add : pd.Series([
                ##########################################50
                str(round(statistics.mean(precision_50)*100,2)) + " ± " +str(round(statistics.stdev(precision_50)*100,2)),
                str(round(statistics.mean(recall_50)*100,2)) + " ± " +str(round(statistics.stdev(recall_50)*100,2)),
                f1_str_50,
                str(round(statistics.mean(ndcg_50)*100,2)) + " ± " +str(round(statistics.stdev(ndcg_50)*100,2)),
                ##########################################n/10
                str(round(statistics.mean(precision_n_10)*100,2)) + " ± " +str(round(statistics.stdev(precision_n_10)*100,2)),
                str(round(statistics.mean(recall_n_10)*100,2)) + " ± " +str(round(statistics.stdev(recall_n_10)*100,2)),
                f1_str_n_10,
                str(round(statistics.mean(ndcg_n_10)*100,2)) + " ± " +str(round(statistics.stdev(ndcg_n_10)*100,2)),
                ##########################################n/4
                str(round(statistics.mean(precision_n_4)*100,2)) + " ± " +str(round(statistics.stdev(precision_n_4)*100,2)),
                str(round(statistics.mean(recall_n_4)*100,2)) + " ± " +str(round(statistics.stdev(recall_n_4)*100,2)),
                f1_str_n_4,
                str(round(statistics.mean(ndcg_n_4)*100,2)) + " ± " +str(round(statistics.stdev(ndcg_n_4)*100,2)),
                ##########################################n/2
                str(round(statistics.mean(precision_n_2)*100,2)) + " ± " +str(round(statistics.stdev(precision_n_2)*100,2)),
                str(round(statistics.mean(recall_n_2)*100,2)) + " ± " +str(round(statistics.stdev(recall_n_2)*100,2)),
                f1_str_n_2,
                str(round(statistics.mean(ndcg_n_2)*100,2)) + " ± " +str(round(statistics.stdev(ndcg_n_2)*100,2)),
                ##########################################n
                str(round(statistics.mean(precision_n)*100,2)) + " ± " +str(round(statistics.stdev(precision_n)*100,2)),
                str(round(statistics.mean(recall_n)*100,2)) + " ± " +str(round(statistics.stdev(recall_n)*100,2)),
                f1_str_n,
                str(round(statistics.mean(ndcg_n)*100,2)) + " ± " +str(round(statistics.stdev(ndcg_n)*100,2))])})

            self.write_pkl_dataset(dataframe=df,name=self.disease_ID)
            return df 

    # =========================== Utilities Functions ========================

    def list_to_pikle(self,list,name):
        """
        Writing list to pickle file
        :param list: input list
        :param: name of the output file
        """
        with open(self.folder_path + 'outputs/' + str(name), 'wb') as f:
            pickle.dump(list, f)

    def read_pickle_list(self,name):
        """
        Reading list from pickle file
        :param: name of the output file
        :return list
        """
        with open(self.folder_path + 'outputs/' + str(name), 'rb') as f:
            tmp_list = pickle.load(f)
        return tmp_list

    def create_empty_dataframe(self,name): 
        """
        Creates an empty dataframe to store the results from the evaluation metrics
        """
        dataframe = {'@': [' ', '  50', ' ', ' ',' ', ' n/10', ' ', ' ',' ', ' n/4', ' ', ' ',' ', ' n/2', ' ', ' ',' ', ' n', ' ', ' ',],
            'Metric': ['P', 'R', 'F1', 'nDCG','P', 'R', 'F1', 'nDCG','P', 'R', 'F1', 'nDCG','P', 'R', 'F1', 'nDCG','P', 'R', 'F1', 'nDCG']
            }  
        dataframe = pd.DataFrame(dataframe)
        dataframe.to_pickle(self.folder_path + 'outputs/pkl_datasets/' + str(name) + ".pkl")

    def read_pkl_dataset(self,name):
        """
        Reading dataset from pickle file
        :param: name of the file
        :return dataframe
        """
        dataframe = pd.read_pickle(self.folder_path + 'outputs/pkl_datasets/' + str(name) + ".pkl")
        return dataframe

    def write_pkl_dataset(self,dataframe,name):
        """
        Writing dataset to pickle file
        :param list: input dataframe
        :param: name of the output file
        """
        self.results_df = dataframe.to_pickle(self.folder_path + 'outputs/pkl_datasets/' + str(name)+".pkl")
        return self.results_df

    def dataframe_to_html(self,dataframe):
        """
        Converting dataframe to html file
        :param dataframe: dataframe -> input to save
        """
        dataframe.to_html(self.folder_path + 'outputs/results_table/' + self.disease_ID +".html",index=True)

    
    def dataframe_to_latex(self,dataframe):
        """
        Encodes input dataframe into latex table
        :param dataframe: dataframe -> input to encode
        """
        print(dataframe.to_latex(index=False)) 