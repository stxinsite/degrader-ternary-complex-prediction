from wepy.hdf5 import WepyHDF5
from wepy.resampling.decisions.clone_merge import MultiCloneMergeDecision
from wepy.boundary_conditions.unbinding import UnbindingBC
from wepy.analysis.contig_tree import ContigTree
from wepy.analysis.parents import ParentForest, resampling_panel
import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
import pickle
import networkx as nx
import sys

sys.setrecursionlimit(9999)
debugnode = (631,27)

job_id = sys.argv[1]

def checknetwork(net):
    for n in net.nodes:
        if 'weights' in net.nodes[n]:
            # look at successors.  
            n_succ = len(list(net.successors(n)))

            #if only one, then the successor should have an equal or greater weight
            if n_succ == 1:
                s = list(net.successors(n))[0]
                if 'weights' in net.nodes[s]:
                    assert net.node[n]['weights'] <= net.node[s]['weights']
            elif n_succ > 1:
                # if more than one, the weight should be equally divided among them
                for s in net.successors(n):
                    if 'weights' in net.nodes[s]:
                        assert net.node[s]['weights'] == net.node[n]['weights']/n_succ
                    
    print("network passed!")

if __name__ == "__main__":

    we = WepyHDF5('6hax_binding_results_REVO.wepy.h5',mode='r')
    
    we.open()

    ncycles = we.num_traj_frames(0, 0)
    nwalkers = we.num_trajs

    print("Building parent forest..")
    contig_tree = ContigTree(we, decision_class=MultiCloneMergeDecision,boundary_condition_class=UnbindingBC)
    contig = contig_tree.span_contig(0)
    pf = ParentForest(contig=contig)

    # set weights for each node
    print("Setting weights..")
    if not os.path.isfile('wt_dict.pkl'):
        wt_dict = {}
        for n in pf.graph.nodes:
            print(n)
            try:
                iter(n)
            except TypeError:
                continue

            if n[0] >= 0 and n[0] < ncycles:
                wt_dict.update({(n[0]-1,n[1]):we.get_traj_field(0,n[1],'weights',frames=[n[0]])})

        with open("wt_dict.pkl","wb") as f:
            pickle.dump(wt_dict,f)
    else:
        with open("wt_dict.pkl","rb") as f:
            wt_dict = pickle.load(f)
    pf.set_node_attributes('weights',wt_dict)
        
    # get the boundary condition data for each node
    print("Setting IRMSD.")
    
    #min_ds = np.array(we.h5['runs/0/progress/min_distances'])
    irmsd = np.loadtxt('csv_files/job_{}_binary_irmsd.csv'.format(job_id), delimiter =',')
    irmsd_dict = {}    
    for n in pf.graph.nodes:
        try:
            iter(n)
        except TypeError:
            continue

        if n[0] >= 0 and n[0] < ncycles:
            irmsd_dict.update({(n[0]-1,n[1]):irmsd[n[0]][n[1]]})


    pf.set_node_attributes('irmsd',irmsd_dict)

    # get a resampling panel
    print("Building resampling panel..")
    res_rec = we.resampling_records([0])
    res_panel = resampling_panel(res_rec)

    # get the warping records
    print("Adding warping records..")
    warp_rec = we.warping_records([0])

    warp_false = {}

    for n in pf.graph.nodes:
        try:
            iter(n)
        except TypeError:
            continue

        if n[0] >= 0 and n[0] < ncycles:
            warp_false.update({(n,False)})
 
    pf.set_node_attributes('warp',warp_false)
    for w in warp_rec:
        pf.graph.node[(w.cycle_idx-1,w.walker_idx)]['warp'] = True

    # initialize exit_wts list
    exit_wts = []
    exit_counts = []
    #cutoffs = np.linspace(0.3,1.0,8)

    # Make list to count number of warping events
    cutoffs = [0.1, 0.15, 0.2]
    n_warp_events = np.zeros_like(cutoffs)

    for c in cutoffs:
        if c == np.min(cutoffs):
            c_idx = 0
        else:
            c_idx += 1
    
        print("Getting exit rate for cutoff = {0}".format(c))

        # reset node weights
        with open("wt_dict.pkl","rb") as f:
            wt_dict = pickle.load(f)
        pf.set_node_attributes('weights',wt_dict)

        # reset node penalties
        penalties = {}
        for n in pf.graph.nodes:
          try:
              iter(n)
          except TypeError:
              continue

          if n[0] >= 0 and n[0] < ncycles:
              penalties.update({(n, False)})


        pf.set_node_attributes('penalty',penalties)

        # shift all nodes forward one cycle?
        mapping = dict([((i,j),(i+1,j)) for i in range(ncycles) for j in range(nwalkers)])
        new_graph = nx.relabel_nodes(pf.graph,mapping)

        binding_counter = 0
        exit_wt = np.zeros((pf.n_steps))
        for cycle in range(ncycles):
            for walker in range(nwalkers):
                node = (cycle,walker)
                if new_graph.has_node((node)):
                    if new_graph.node[node]['irmsd'] < c:
                        #print(new_graph.node[node]['penalty'])
                        w = new_graph.node[node]['weights'][0][0] - new_graph.node[node]['penalty']
                        #print(w)
                        assert w >= 0
                            #print("exit point with weight {0}, penalty {1}".format(new_graph.node[node]['weights'][0][0],new_graph.node[node]['penalty']))
                        if w > 0:
                            binding_counter += 1
                        exit_wt[node[0]] += w
                        new_graph.node[node]['penalty'] = w

                    # if you warped, your penalty disappears
                    if new_graph.node[node]['warp'] == True:
                        n_warp_events[c_idx] += 1
                        new_graph.node[node]['penalty'] = 0
                    
                    if node[0]+1 < ncycles:
                        # if you were squashed, subtract your penalty from the squasher
                        if res_panel[node[0]][0][node[1]][0] == 3:
                            # squashed!
                            squashed_into = (node[0]+1, res_panel[node[0]][0][node[1]][1][0])
                            if squashed_into == debugnode:
                                print("Adding weight ({0}) to penalty for node {1} (weight = {2}): new penalty: {3}".format(
                                    new_graph.node[node]['penalty'], squashed_into, new_graph.node[squashed_into]['weights'],
                                    new_graph.node[squashed_into]['penalty']+new_graph.node[node]['penalty']))
                                print("sq from {0}".format(node))
                            new_graph.node[squashed_into]['penalty'] += new_graph.node[node]['penalty']
                            assert new_graph.node[squashed_into]['penalty'] <= new_graph.node[squashed_into]['weights']

                        # otherwise, pass the penalty on to your successors
                        n_succ = len(list(new_graph.successors(node)))
                        for s in new_graph.successors(node):                            
                            if s == debugnode:
                                print("Adding weight ({0}) to penalty for node {1} (weight = {2}): new penalty: {3}".format(
                                    new_graph.node[node]['penalty']/n_succ, s, new_graph.node[s]['weights'],
                                    new_graph.node[s]['penalty']+new_graph.node[node]['penalty']/n_succ))
                                print("succession from {0}".format(node))
                            new_graph.node[s]['penalty'] += new_graph.node[node]['penalty']/n_succ
                            assert new_graph.node[s]['penalty'] <= new_graph.node[s]['weights']

        # save exit_wt to array
        exit_wts.append(exit_wt)
        exit_counts.append([binding_counter])


    exit_wt_totals = [np.zeros(shape=e.shape) for e in exit_wts]
    for i,e in enumerate(exit_wts):
        for j,e_val in enumerate(e):
            exit_wt_totals[i][j:] += e_val

            
    pickle.dump(exit_wt_totals, open('job_{}_exit_weights.pkl'.format(job_id), 'wb'))
    pickle.dump(exit_counts, open('job_{}_exit_counts.pkl'.format(job_id), 'wb'))                                                                                                                      
    we.close()
