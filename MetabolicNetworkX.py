import networkx as nx
import xml.etree.ElementTree as ET
from pyvis.network import Network
import matplotlib.pyplot as plt


class MetabolicNetworkX():

    def __init__(self, split_rev=True):
        self.mnetwork = nx.DiGraph()  # networkX data structure used to implement Metabolic Networks
        self.split_rev = split_rev

    def print_network(self):
        '''
        simple node print statement

        :return: print'node -> sucessors'
        '''
        for node in self.mnetwork.nodes:
            print(node, '->', list(self.mnetwork.successors(node)))

    def loadXML(self, filepath):
        '''
        loads metabolic network from XML-SBML file;
        tested with microrganism files from BioModels - EMBL-EBI.
        TO DO: clean method, implement formal SBML reader

        :param filepath: (str) path to XML-SBML file
        :return: overrides self.mnetwork with metabolic network structure
        '''
        datatree = ET.parse(filepath)
        root = datatree.getroot()
        rr = root[0].findall('{http://www.sbml.org/sbml/level3/version1/core}listOfReactions')
        for reaction in rr[0]:
            child = reaction.getchildren()
            try:
                if len(child) > 1:

                    listOfReactants = child[[child.index(i) for i in child if '}listOfReactants' in str(i)][0]]
                    for reagent in listOfReactants.getchildren():

                        if reaction.attrib['reversible'] == 'true':
                            self.mnetwork.add_node(reagent.attrib['species'], data="metabolite")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_f', data="reaction")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_b', data="reaction")

                            self.mnetwork.add_edge(reagent.attrib['species'], reaction.attrib['id'] + '_f')
                            self.mnetwork.add_edge(reaction.attrib['id'] + '_b', reagent.attrib['species'])

                        else:
                            self.mnetwork.add_node(reagent.attrib['species'], data="metabolite")
                            self.mnetwork.add_node(reaction.attrib['id'], data="reaction")

                            self.mnetwork.add_edge(reagent.attrib['species'], reaction.attrib['id'])

                    listOfProducts = child[[child.index(i) for i in child if '}listOfProducts' in str(i)][0]]
                    for product in listOfProducts.getchildren():

                        if reaction.attrib['reversible'] == 'true':

                            self.mnetwork.add_node(product.attrib['species'], data="metabolite")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_f', data="reaction")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_b', data="reaction")

                            self.mnetwork.add_edge(reaction.attrib['id'] + '_f', product.attrib['species'])
                            self.mnetwork.add_edge(product.attrib['species'], reaction.attrib['id'] + '_b')
                        else:

                            self.mnetwork.add_node(product.attrib['species'], data="metabolite")
                            self.mnetwork.add_node(reaction.attrib['id'], data="reaction")

                            self.mnetwork.add_edge(reaction.attrib['id'], product.attrib['species'])



                else:
                    if reaction.attrib['reversible'] == 'true':
                        for M in child[0]:
                            self.mnetwork.add_node(M.attrib['species'], data="metabolite")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_f', data="reaction")
                            self.mnetwork.add_node(reaction.attrib['id'] + '_b', data="reaction")

                            self.mnetwork.add_edge(M.attrib['species'], reaction.attrib['id'] + '_f')
                            self.mnetwork.add_edge(reaction.attrib['id'] + '_b', M.attrib['species'])

            except:
                print(reaction.attrib['id'], 'couln\'t be imported to graph')
                print(reaction.getchildren())
                continue

    def load_from_file(self, filepath):
        '''
        loads metabolic network from TXT file;
        file should use the following format:
        REACTION: METABOLITE => PRODUCT
        note that reversible reactions are represented by '<=>',
        and multiple metabolites and prodcts are seperated by ' + '

        :param filepath: (str) path to XML-SBML file
        :return: verrides self.mnetwork with metabolic network structure
        '''
        rf = open(filepath)
        for line in rf:
            if ":" in line:
                tokens = line.split(":")
                reac_id = tokens[0].strip()
                self.mnetwork.add_node(reac_id, data="reaction")
                rline = tokens[1]
            else:
                raise Exception("Invalid line:")
            if "<=>" in rline:
                left, right = rline.split("<=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in self.mnetwork.graph:
                        self.mnetwork.add_node(met_id, data="metabolite")
                    if self.split_rev:
                        self.mnetwork.add_node(reac_id + "_b", data="reaction")
                        self.mnetwork.add_edge(met_id, reac_id)
                        self.mnetwork.add_edge(reac_id + "_b", met_id)
                    else:
                        self.mnetwork.add_edge(met_id, reac_id)
                        self.mnetwork.add_edge(reac_id, met_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in self.mnetwork.graph:
                        self.mnetwork.add_node(met_id, data="metabolite")
                    if self.split_rev:
                        self.mnetwork.add_edge(met_id, reac_id + "_b")
                        self.mnetwork.add_edge(reac_id, met_id)
                    else:
                        self.mnetwork.add_edge(met_id, reac_id)
                        self.mnetwork.add_edge(reac_id, met_id)
            elif "=>" in line:
                left, right = rline.split("=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in self.mnetwork.graph:
                        self.mnetwork.add_node(met_id, data="metabolite")
                    self.mnetwork.add_edge(met_id, reac_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in self.mnetwork.graph:
                        self.mnetwork.add_node(met_id, data="metabolite")
                    self.mnetwork.add_edge(reac_id, met_id)
            else:
                raise Exception("Invalid line:")

    def summary(self):
        '''
        Provides basic information on the network,
        Number of nodes and number of edges

        :return: prints newline seperated strings
        '''
        print('Network has:', self.mnetwork.number_of_nodes(), 'nodes')
        print('Network has:', self.mnetwork.number_of_edges(), 'edges\n')

    def see_disconnected(self):
        '''
        importing networks leads to possible data fragments and isolated nodes or graphs
        this method allows to see all subgraphs and isolated nodes in the data structure

        :return: prints all disconneted graphs within mnetwork
        '''
        for subgraph in list(nx.weakly_connected_components(self.mnetwork)):
            print(subgraph)

    def remove_disconnected(self):
        '''
        importing networks leads to possible data fragments and isolated nodes or graphs,
        these fragments can be due to errors or leftovers in the original file or bugs on import
        some class methods can't equate these fragments and metrics can perform diferently in their presence
        this method removes all graphs and nodes isolated from the main branch

        :return: overides self.mnetwork removing disconected nodes and subgraphs
        '''
        for subgraph in list(nx.weakly_connected_components(self.mnetwork)):
            listed_nodes = list(subgraph)
            if len(listed_nodes) < (self.mnetwork.number_of_edges() * 0.1):
                self.mnetwork.remove_nodes_from(listed_nodes)

    def see_graph(self, htmlname='MetabolicNetwork.html', height="1080px", width="1080px", colorHighNodes=5):
        '''
        #***DEPRECATED METHOD***#
        #***use see_graph2***#
        Provides a visual representation of the network using pyvis
        Metabolites are represented by yellow dots
        Reationss are represented by green triangles
        Metabolite Highlights are represented by red dots and denote highest in_degrees metabolites
        if nodes are not metabolites or reaction default blue dot is used

        :param htmlname:(str) name of exit file
        :param height:(num) height of html window display
        :param width:(num) width of html window display
        :param colorHighNodes:(int) number of nodes to color highlight, if 0 will not highlight any node

        :return: creates html file and opens it on default browser
        '''
        nt = Network(height, width)
        nt.from_nx(self.mnetwork)

        for dit in nt.nodes:
            if str(dit['title']).startswith('R'):
                dit['shape'] = 'triangle'
                dit['color'] = 'green'
            if str(dit['title']).startswith('M'):
                dit['color'] = 'rgb(255, 211, 0)'
        nt.show_buttons(filter_=['physics', 'edges'])
        nt.inherit_edge_colors_from(False)

        if colorHighNodes > 0:
            mustcolor = self.biggestM_indegrees(colorHighNodes)
            for dit2 in nt.nodes:
                if str(dit2['title']) in mustcolor:
                    dit2['color'] = 'red'

        nt.barnes_hut(gravity=-3000)
        nt.toggle_stabilization(True)

        nt.show(htmlname)

    def see_graph2(self, htmlname='MetabolicNetwork.html', highlightCentrality=False, height="1080px", width="1080px",
                   h_degree=True, cent_between=True, cent_close=True, top_num=10):
        '''
        Provides a visual representation of the network using pyvis,
            allowing topologic metric highlights on top_num nodes
        Metabolites are represented by yellow dots
        Reations are represented by green triangles
        if nodes are not metabolites or reaction default blue dot is used


        :param htmlname:(str) name of exit file
        :param highlightCentrality: (bool) When set to True highlights nodes and eges based on topological metrics.
                                Default is set to False
        :param height: (num) height of html window display
        :param width: (num) width of html window display
        :param h_degree: (bool) Isn't used unless highlightCentrality is True.
                                Computes Highest Degree Metric for highlight.
                                Default is set to True
        :param cent_between: (bool) Isn't used unless highlightCentrality is True.
                                Computes Betweenness Centrality for highlight.
                                Default is set to True
        :param cent_close: (bool) Isn't used unless highlightCentrality is True.
                                Computes Closeness Centrality for highlight.
                                Default is set to True
        :param top_num: (int) Number of nodes to highlight for each selected metric. Default is set to 10

        :return: creates html file and opens it on default browser
        '''

        nt = Network(height, width)
        nt.from_nx(self.mnetwork)
        nt.options.configure = True
        nt.options.edges.color = 'rgb(200,200,200)'  # lightgray-ish

        if highlightCentrality:
            print('You have selected to highlight nodes and edges based on their centrality measures\n')

        # basic colour applier
        for dit in nt.nodes:
            if str(dit['title']).startswith('R'):
                dit['shape'] = 'triangle'
                dit['color'] = 'green'
            if str(dit['title']).startswith('M'):
                dit['color'] = 'rgb(255, 211, 0)'  # yellow-ish

        if highlightCentrality and h_degree:  # apply color to highest degree nodes
            to_col = self.highest_degrees(top=top_num)
            tuple_col = sghelp_listMR(to_col)
            print('Highlighted Nodes based on highest degree:\n', tuple_col)
            to_col_M = tuple_col[0]
            to_col_R = tuple_col[1]

            for tm in to_col_M:
                for node1 in nt.nodes:
                    if str(node1['title']) == tm:
                        node1['color'] = 'red'

            for tr in to_col_R:
                for node2 in nt.nodes:
                    if str(node2['title']) == tr:
                        node2['color'] = 'rgb(0,0,139)'  # darkblue-ish

        if highlightCentrality and cent_close:  # apply size change to highest closseness nodes
            to_fat = self.highest_closeness(top=top_num)
            print('Highlighted Nodes based on Closeness Centrality:\n', to_fat)
            sizevalue = 12 * top_num

            for tf in to_fat:
                for node3 in nt.nodes:
                    if str(node3['title']) == tf:
                        node3['size'] = sizevalue
                        sizevalue -= top_num

        if highlightCentrality and cent_between:  # apply width change to edges of nodes with highest betweeness
            to_thick = self.highest_betweenness(top=top_num)
            print('Highlighted Edges based on Betweenness Centrality:\n', to_thick)

            for tt in to_thick:
                for edge in nt.edges:
                    if edge['from'] == tt:
                        edge['value'] = 1.5
                        for node in nt.nodes:
                            if node['title'] == edge['from']:
                                edge['color'] = node['color']

        nt.barnes_hut(gravity=-3000)
        nt.toggle_stabilization(True)
        nt.show_buttons(filter_=['physics'])
        nt.show(htmlname)

    def biggestM_indegrees(self, n=5):
        '''
        #***DEPRECATED METHOD***#
        #***use in_degree_list, all_degree_list or all_degreetype_dict***#
        calculates highest in_deree nodes
        :param n: (int) number of nodes in list

        :return: (list) sorted metabolits with highest in_degree
        '''
        allindeg = list(self.mnetwork.in_degree())
        res = []
        ll = sorted(allindeg, key=lambda x: x[1])
        for i in range(len(ll) - 1, len(ll) - n - 1, -1):
            if str(ll[i][0]).startswith('M'): res.append(str(ll[i][0]))
        return res

    def remove_reactions(self):
        '''
        Default import creates a Metabolite-Reaction Bipartite graph
        this method removes Reaction and Creates a Metabolite-Metabolite graph
        TO DO: revert change and implement Reaction equivalent

        :return: overrides self.mnetwork
        '''
        rem = []
        for node1 in self.mnetwork.nodes:
            if (self.mnetwork.node[node1]['data'] == 'metabolite'):
                sucs = list(self.mnetwork.successors(node1))
                for s in sucs:
                    sucs_r = self.mnetwork.successors(s)
                    for s2 in sucs_r:
                        if node1 != s2:
                            self.mnetwork.add_edge(node1, s2)
            else:
                rem.append(node1)
        self.mnetwork.remove_nodes_from(rem)

    # ---- Degree Distribution
    def in_degree_list(self):
        '''
        calculates in degrees for all nodes
        :return: (list) [(Node,in_degree value)] sorted
        '''
        return sorted(self.mnetwork.in_degree(), key=lambda x: x[1])

    def out_degree_list(self):
        '''
        calculates out degrees for all nodes
        :return: (list) [(Node,out_degree value)] sorted
        '''
        return sorted(self.mnetwork.out_degree(), key=lambda x: x[1])

    def all_degree_list(self):
        '''
        calculates degrees for all nodes
        :return: (list) [(Node,degree value)] sorted
        '''
        return sorted(self.mnetwork.degree(), key=lambda x: x[1])

    def all_degreetype_dict(self, deg_type='inout'):
        '''
        calculates (in/out)degrees for all nodes
        :param deg_type: (str) degree type. Default is inout (all degrees)

        :return: (dict) {node: degree_value}
        '''

        if deg_type == 'inout':
            dg = self.mnetwork.degree()

        elif deg_type == 'in':
            dg = self.mnetwork.in_degree()

        else:
            dg = self.mnetwork.out_degree()

        return tuplelist_to_dict(dg)

    def mean_degree(self, deg_type="inout"):
        '''
        Calculates mean degree
        :param deg_type: (str) degree type. Default is inout (all degrees)

        :return: (num) Mean degree value
        '''
        degs = self.all_degreetype_dict(deg_type)
        return sum(degs.values()) / float(len(degs))

    def prob_degree(self, deg_type="inout"):
        '''
        Calculates probability of degree value in the network
        :param deg_type: (str) degree type. Default is inout (all degrees)

        :return: (dict) {degree: probability}
        '''
        degs = self.all_degreetype_dict(deg_type)
        res = {}
        for k in degs.keys():
            if degs[k] in res.keys():
                res[degs[k]] += 1
            else:
                res[degs[k]] = 1
        for k in res.keys():
            res[k] /= float(len(degs))
        return res

    # ---- Clustering
    def mean_clustering_perdegree(self, deg_type="inout"):
        '''
        calculates clustering coefecient for all degrees (C(k))
        :param deg_type: (str) degree type. Default is inout (all degrees)

        :return: (dict) {degree: cluestering_coef}
        '''
        degs = self.all_degreetype_dict(deg_type)
        ccs = nx.clustering(self.mnetwork)
        degs_k = {}
        for k in degs.keys():
            if degs[k] in degs_k.keys():
                degs_k[degs[k]].append(k)
            else:
                degs_k[degs[k]] = [k]
        ck = {}
        for k in degs_k.keys():
            tot = 0
            for v in degs_k[k]: tot += ccs[v]
            ck[k] = float(tot) / len(degs_k[k])
        return ck

    # ---- Centrality Measures
    def highest_degrees(self, top=10):
        '''
        Calculates nodes with highest degrees in network
        :param top: (int) number of nodes to display

        :return: (list) sorted high->low Top Nodes with highest degree
        '''
        all_deg = self.all_degree_list()
        ord_deg = sorted(all_deg, key=lambda x: x[1], reverse=True)
        return list(map(lambda x: x[0], ord_deg[:top]))

    def betweenness_centrality(self):
        '''
        Calculates betweenness centrality metric for each node
        :return: (dict) {node:value}
        '''
        return nx.algorithms.betweenness_centrality(self.mnetwork)

    def highest_betweenness(self, top=10):
        '''
        Calculates nodes with highest betweenness in network
        :param top: (int) number of nodes to display

        :return: (list) sorted high->low Top Nodes with highest betweenness
        '''
        bb = self.betweenness_centrality()
        ord_bb = sorted(list(bb.items()), key=lambda x: x[1], reverse=True)
        return list(map(lambda x: x[0], ord_bb[:top]))

    def closeness_centrality(self, nd=None):
        '''
        Calculates closeness centrality metric for each node
        :param nd: (str) if used calculates metric for specific node. Default is None

        :return: (dict) {node:value}
        '''
        return nx.algorithms.closeness_centrality(self.mnetwork, u=nd)

    def highest_closeness(self, top=10):
        '''
        Calculates nodes with highest closeness in network
        :param top: (int) number of nodes to display

        :return: (list) sorted high->low Top Nodes with highest closeness
        '''
        cc = self.closeness_centrality()
        ord_cl = sorted(list(cc.items()), key=lambda x: x[1], reverse=True)
        return list(map(lambda x: x[0], ord_cl[:top]))

    # ----
    def topological_analysis_full(self):
        '''
         Calculates and displays all implemented network analysis
        :return: prints topological data analysis
        '''
        print('###Network Topological Analysis###\n')
        # -------
        print('1. Degree Distribution\n')
        print('\tNetwork Mean Degree:')
        print('\t\t Mean InDegree:', self.mean_degree('in'))
        print('\t\t Mean OutDegree:', self.mean_degree('out'))
        print('\t\t Mean (Inout)Degree:', self.mean_degree())
        print()
        print('\tNetwork Probability (Inout)Degree:\n', self.prob_degree())
        self.plot_prob_degree(title='Network Probability Degree', xlabel='k', ylabel='P(k)')
        print()
        print()
        # -------
        print('2. Shortest Path Analysis')
        print('\tMean Distances:', nx.average_shortest_path_length(self.mnetwork))
        print()
        print()
        # -------
        print('3. Clustering Coefficients')
        print('\tAverage Clustering:', nx.average_clustering(self.mnetwork))
        print('\tMean Clustering per Degree:\n', self.mean_clustering_perdegree())
        self.plot_clustering_perdegree(title='Network Clustering Per Degree', xlabel='k', ylabel='C(k)')
        print()
        print()
        # -------
        print('4. Hubs and Centrality Measures')
        print('\tNodes w/ Highest Degrees:', self.highest_degrees())
        print('\tNodes w/ Highest Betweenness:', self.highest_betweenness())
        print('\tNodes w/ Highest Closeness:', self.highest_closeness())

    def plot_prob_degree(self, logDegrees=False, logProb=False, title='', xlabel='', ylabel='', color='g',
                         xlim=(None, None),
                         ylim=(None, None)):
        '''
        Plots probability degree of network
        :param logDegrees: (bool) If active logarithmizes Degree Values. Default is False
        :param logProb: (bool) If active logarithmizes Probability Values. Default is False
        :param title: (str) Title of the plot
        :param color: (str) Matplot Lib color format
        :param xlim: (tuple)  Set the x limits of the axes. Default is no limits
        :param ylim: (tuple) Set the y limits of the axes. Default is no limits

        :return: matplotlib plot
        '''
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        pd = self.prob_degree()
        x = list(pd.keys())
        y = list(pd.values())

        if logDegrees:
            x = log_list_num(x)

        if logProb:
            y = log_list_num(y)

        plt.plot(x, y, 'ro', color=color)

        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        plt.title(title)
        plt.show()

    def plot_clustering_perdegree(self, logDegrees=False, logClust=False, title='', xlabel='', ylabel='', color='g',
                                  xlim=(None, None),
                                  ylim=(None, None)):
        '''
        plots cluestering coefecient per degree
        :param logDegrees: (bool) If active logarithmizes Degree Values. Default is False
        :param logClust: (bool) If active logarithmizes Cluster Coef Values. Default is False
        :param title: (str) Title of the plot
        :param color: (str) Matplot Lib color format
        :param xlim: (tuple)  Set the x limits of the axes. Default is no limits
        :param ylim: (tuple) Set the y limits of the axes. Default is no limits

        :return: matplotlib plot
        '''
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        cl = self.mean_clustering_perdegree()
        x = list(cl.keys())
        y = list(cl.values())

        if logDegrees:
            x = log_list_num(x)

        if logClust:
            y = log_list_num(y)

        plt.plot(x, y, 'ro', color=color)

        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        plt.title(title)
        plt.show()


########################################################################################################################


def tuplelist_to_dict(tuplelist):
    res = {}
    for node, value in tuplelist:
        res[node] = value
    return res


def is_in_tuple_list(tl, val):
    res = False
    for (x, y) in tl:
        if val == x: return True
    return res


def sghelp_listMR(nodelist):  # see_graph2 helper, seperates Metabolites from Reactions
    metabolits = []
    reactions = []
    for node in nodelist:
        if str(node).startswith('M'):
            metabolits.append(node)
        else:
            reactions.append(node)

    return (metabolits, reactions)


def log_list_num(lista):
    from math import log10
    return [log10(y) for y in lista]


def make_probdegree_plot(graph, title='', color='g'):
    pd = graph.prob_degree()
    # print(pd)
    # pdo = sorted( list(pd.items()), key=lambda x : x[0])
    plt.plot(list(pd.keys()), list(pd.values()), 'ro', color=color)
    # plt.xlim(0, 50)
    # plt.ylim(0,0.01)
    plt.title(title)
    plt.show()


def make_clustering_plot(graph, title='', color='g'):
    cl = graph.mean_clustering_perdegree()
    plt.plot(list(cl.keys()), list(cl.values()), 'ro', color=color)
    plt.title(title)
    plt.show()


def testXML():
    T = MetabolicNetworkX()
    T.loadXML('MODEL1507180063_url.xml')
    T.see_graph('modeloxml.html')


def testTXT():
    T = MetabolicNetworkX()
    T.load_from_file('ecoli.txt')
    T.see_graph('modelotxt.html')


def testMetG():
    T = MetabolicNetworkX()
    T.load_from_file('ecoli.txt')

    print(T.mnetwork.size())

    T.remove_reactions()
    T.see_graph('sometabolitos.html')


def testDegreeDist():
    T = MetabolicNetworkX()
    T.load_from_file('ecoli.txt')
    # print(T.in_degree_list())
    # print(T.out_degree_list())
    # print(T.all_degree_list())
    # print(T.all_degreetype_dict())
    print(T.mean_degree())

    pd = T.prob_degree()
    # pdo = sorted( list(pd.items()), key=lambda x : x[0])
    plt.bar(list(pd.keys()), pd.values(), color='g')
    plt.xlim(0, 50)
    plt.show()


def testDegreeComp():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')

    ST = MetabolicNetworkX()
    ST.loadXML('MODEL1507180063_url.xml')  # Streptococcus thermophilus
    ST.remove_disconnected()

    SF = MetabolicNetworkX()
    SF.mnetwork = nx.scale_free_graph(EC.mnetwork.number_of_nodes(), seed=2424)

    R = MetabolicNetworkX()
    R.mnetwork = nx.gnp_random_graph(EC.mnetwork.number_of_nodes(), 0.5, seed=2424)  # Erdös–Rényi (ER) model

    print('Comparing Theoretical Scale free and Random  Networks with Real Biological Data')

    print('Mean Degree')
    print(' E.coli: ', EC.mean_degree(), '\n',
          'S.thermophilus: ', ST.mean_degree(), '\n',
          'Theoretical Scale Free:', SF.mean_degree(), '\n',
          'Theoretical Random:', R.mean_degree())

    make_probdegree_plot(EC, 'Ecoli')
    make_probdegree_plot(ST, 'Stermo', 'r')
    make_probdegree_plot(SF, 'Theoretical Scale Free', 'b')
    make_probdegree_plot(R, 'Theoretical Random', 'm')

    # foi corrido uma vez apenas para construir html
    # G.see_graph('Ecoli.html')
    # SF.see_graph('Scalefree.html',colorHighNodes=0)
    # R.see_graph('RandomGen.html',colorHighNodes=0)


def testShortPath():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')

    ST = MetabolicNetworkX()
    ST.loadXML('MODEL1507180063_url.xml')
    ST.remove_disconnected()  # limpar residuos

    SF = MetabolicNetworkX()
    SF.mnetwork = nx.scale_free_graph(EC.mnetwork.number_of_nodes(), seed=2424)

    R = MetabolicNetworkX()
    R.mnetwork = nx.gnp_random_graph(EC.mnetwork.number_of_nodes(), 0.5, seed=2424)

    print('###############')
    print(' E.coli: ', nx.average_shortest_path_length(EC.mnetwork), '\n',
          'S.thermophilus: ', nx.average_shortest_path_length(ST.mnetwork), '\n',
          'Theoretical Scale Free:', nx.average_shortest_path_length(SF.mnetwork), '\n',
          'Theoretical Random:', nx.average_shortest_path_length(R.mnetwork))


def testClustering():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')

    ST = MetabolicNetworkX()
    ST.loadXML('MODEL1507180063_url.xml')
    ST.remove_disconnected()  # limpar residuos

    SF = MetabolicNetworkX()
    SF.mnetwork = nx.scale_free_graph(EC.mnetwork.number_of_nodes(), seed=2424)

    R = MetabolicNetworkX()
    R.mnetwork = nx.gnp_random_graph(EC.mnetwork.number_of_nodes(), 0.5, seed=1111)

    # print(nx.average_clustering(EC.mnetwork))
    # print(nx.average_clustering(ST.mnetwork))
    # print(nx.average_clustering(SF.mnetwork))
    # print(nx.average_clustering(R.mnetwork))

    print('###############')
    print('Clustering Analysis')
    print(' Average Clustering E.coli: ', nx.average_clustering(EC.mnetwork), '\n',
          'Average Clustering S.thermophilus: ', nx.average_clustering(ST.mnetwork), '\n',
          # 'Average Clustering Theoretical Scale Free:', nx.average_clustering(SF.mnetwork), '\n',
          'Average Clustering Theoretical Random:', nx.average_clustering(R.mnetwork))

    make_clustering_plot(EC, 'Ecoli')
    make_clustering_plot(ST, 'Stermo', 'r')
    # make_clustering_plot(SF, 'Theoretical Scale Free', 'b')
    make_clustering_plot(R, 'Theoretical Random', 'm')


def testCentrality():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')
    print('Highest Degree:\n', EC.highest_degrees())
    print()
    print('Betweenness Centrality:\n', EC.betweenness_centrality())
    print()
    print('Closeness Centrality:\n', EC.closeness_centrality())


def testCentralityVisual():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')
    EC.see_graph2('testCentrality.html', highlightCentrality=True)


def testCentralityVisualSpec():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')
    EC.see_graph2('AAtestCentrality.html',
                  highlightCentrality=True,
                  h_degree=True,
                  cent_between=False,
                  cent_close=True,
                  top_num=10)


def testFullTopologicalAnalysis():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')
    EC.topological_analysis_full()


def testGraphs():
    EC = MetabolicNetworkX()
    EC.load_from_file('ecoli.txt')
    EC.plot_prob_degree(title='E.coli Prob Degree', logDegrees=True, ylim=(0, 0.01))
    EC.plot_clustering_perdegree(color='r')

def makenetworks():
    AN = MetabolicNetworkX() #Aspergillus niger
    AN.loadXML('MODEL1507180047_url.xml')
    AN.remove_disconnected()

    ST = MetabolicNetworkX() #Streptococcus thermophilus
    ST.loadXML('MODEL1507180063_url.xml')
    ST.remove_disconnected()

    SA = MetabolicNetworkX() #Staphylococcus aureus
    SA.loadXML('MODEL1507180070_url.xml')
    SA.remove_disconnected()

    AN.see_graph2('Aniger.html')
    AN.see_graph2('Aniger_withcentrality.html',highlightCentrality=True)

    ST.see_graph2('Stermo.html')
    ST.see_graph2('Stermo_withcentrality.html',highlightCentrality=True)

    SA.see_graph2('Saureus.html')
    SA.see_graph2('Saureus_withcentrality.html',highlightCentrality=True)

# testXML()
# testTXT()
# testMetG()
# testDegreeDist()
# testDegreeComp()
# testShortPath()
# testClustering()
# testCentrality()
# testCentralityVisual()
# testCentralityVisualSpec()
# testFullTopologicalAnalysis()
# testGraphs()
# makenetworks()