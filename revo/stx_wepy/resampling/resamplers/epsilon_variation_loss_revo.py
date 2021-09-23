import random as rand
import itertools as it
import logging
import numpy as np

# TODO: rename to PatchedREVOResampler and refactor to only necessary changes
from stx_wepy.resampling.resamplers.revo import REVOResampler

class EpsilonVariationLossREVOResampler(REVOResampler):
    BEST_PROG_METHODS = ('min', 'max')

    def __init__(self,
                 bc_condition=None,
                 epsilon=np.inf,
                 best_prog =None,
                 **kwargs):


        """Constructor for the REVO Resampler.

        Parameters
        ----------

        See REVO Resampler parameters

        bc_condition: wepy boundary condition class.
            The boundary condition to used to calculate simulation
            progress.

        epsilon: float
            The range along the progress coordinate to clone walkers.
            Walkers eligible for cloning is is the max progress walker + epsilon. 

        best_prog : str (either min or max)
            Determining if the progress is maximizing or minimizing the
            progress coordinate.

        """


        super().__init__(**kwargs)

        assert bc_condition is not None, "Boundary Condition must be given."
        assert best_prog in self.BEST_PROG_METHODS, 'best_prog is not valid. Valid options are: {}'.format(self.BEST_PROG_METHODS)
        assert epsilon >= 0, 'epsilon must be a positive value.'

        self.bc_cond = bc_condition

        self.epsilon = epsilon

        self.best_prog = best_prog
                
    def _calc_variation_loss(self, walker_variation, weights, eligible_pairs):
        """Calculates the loss to variation through merging of eligible walkers.                                                                   

        Parameters                                                                                                                                      
        ----------                                                                                                                       

        walker_variations : arraylike of shape (num_walkers)                                                                                            
           The Vi value of each walker.                                                                                                                       

        weights : list of float                                  
            The weights of all walkers. The sum of all weights should be 1.0.

        eligible_pairs : list of tuples
            Pairs of walker indexes that meet the criteria for merging.
                                                                                                                                                          
        Returns                                                                                                                              
        -------                                                       

        variation_loss_list : tuple
            A tuple of the walker merge pair indicies that meet the criteria 
            for merging and minimize variation loss.
         """

        v_loss_min = np.inf
        
        min_loss_pair = ()

        for pair in eligible_pairs:
            walker_i = pair[0]
            walker_j = pair[1]

            wt_i = weights[walker_i]
            wt_j = weights[walker_j]

            v_i = walker_variation[walker_i]
            v_j = walker_variation[walker_j]

            v_loss = (wt_j * v_i + wt_i * v_j) / (wt_i + wt_j)

            if v_loss < v_loss_min:
                min_loss_pair = pair

                v_loss_min = v_loss

        return min_loss_pair

    def _find_eligible_merge_pairs(self, weights, distance_matrix, max_var_idx, num_walker_copies):
        """ Find pairs of walkers that are eligible to be merged.

        Parameters
        ----------

        weights : list of float                                                                                                 
            The weights of all walkers. The sum of all weights should be 1.0.

        distance_matrix : list of arraylike of shape (num_walkers)
            The distance between every walker according to the distance metric.

        max_var_idx : float 
            The index of the walker that had the highest walker variance
            and is a candidate for cloning.

        num_walker_copies : list of int                                                                                                              
            The number of copies of each walker.                                                                                                       
             0 means the walker is not exists anymore.                                                                                                  
             1 means there is one of the this walker.                                                                                                 
             >1 means it should be cloned to this number of walkers.                                                                                        

        Returns
        -------
        
        eligible_pairs : list of tuples                                                                                                                             
            Pairs of walker indexes that meet the criteria for merging.

        """

        eligible_pairs = []

        for i in range(len(weights) - 1):
            for j in range(i + 1, len(weights)):
                if i != max_var_idx and j != max_var_idx:
                    if num_walker_copies[i] == 1 and num_walker_copies[j] == 1:
                        if weights[i] + weights[j] < self.pmax:
                            if distance_matrix[i][j] < self.merge_dist:
                                eligible_pairs.append((i,j))

        return eligible_pairs

    def _calc_progress(self, walkers):

        """ Determine the progress of each walker toward the boundary condition.

        Parameters
        ----------

        walkers: list of walkers

        Returns
        ------
        
        progress: list of floats
            list containing the progress of each walker.
        """

        n_walkers = len(walkers)

        progress = np.zeros([n_walkers])
        
        for walker in range(n_walkers):
            warp, prog_dic = self.bc_cond._progress(walkers[walker])
            keys = [*prog_dic][0]
            progress[walker] = prog_dic[keys]


        return progress

            

    def decide(self, walker_weights, num_walker_copies, distance_matrix, walker_progress):
        """Optimize the trajectory variation by making decisions for resampling.

        Parameters
        ----------

        walker_weights : list of flaot
            The weights of all walkers. The sum of all weights should be 1.0.

        num_walker_copies : list of int
            The number of copies of each walker.
            0 means the walker is not exists anymore.
            1 means there is one of the this walker.
            >1 means it should be cloned to this number of walkers.

        distance_matrix : list of arraylike of shape (num_walkers)

        walker_progress: list of floats of len num_walkers)
            A list of walker progress toward the boundary condition.

        Returns
        -------
        variation : float
            The optimized value of the trajectory variation.

        resampling_data : list of dict of str: value
            The resampling records resulting from the decisions.

        """

        if self.best_prog == 'min':
            threshold = np.min(walker_progress) + self.epsilon

        elif self.best_prog == 'max':
            threshold = np.max(walker_progress) - self.epsilon

        else:
            raise Exception ('best_prog needs to be string containing either "min" or "max".')
        
        num_walkers = len(walker_weights)

        variations = []
        merge_groups = [[] for i in range(num_walkers)]
        walker_clone_nums = [0 for i in range(num_walkers)]

        # make copy of walkers properties
        new_walker_weights = walker_weights.copy()
        new_num_walker_copies = num_walker_copies.copy()


        # calculate the initial variation which will be optimized
        variation, walker_variations = self._calcvariation(walker_weights,
                                                           new_num_walker_copies,
                                                           distance_matrix)
        variations.append(variation)

        # maximize the variance through cloning and merging
        logging.info(f"Starting variance optimization: {variation}")

        productive = True
        while productive:
            productive = False
            # find min and max walker_variationss, alter new_amp

            # initialize to None, we may not find one of each
            min_idx = None
            max_idx = None

            # selects a walker with minimum walker_variations and a walker with
            # maximum walker_variations walker (distance to other walkers) will be
            # tagged for cloning (stored in maxwind), except if it is
            # already a keep merge target
            max_tups = []
            for i, value in enumerate(walker_variations):
                # 1. must have an amp >=1 which gives the number of clones to be made of it
                # 2. clones for the given amplitude must not be smaller than the minimum probability
                # 3. must not already be a keep merge target
                # 4. Must be within the progress coordinate threshold

                # Did the walker get merged away in a previous resampling step
                if new_num_walker_copies[i] >= 1:

                    # Is the walker weight after cloning greater than pmin
                    if new_walker_weights[i]/(new_num_walker_copies[i] + 1) > self.pmin:
                        if len(merge_groups[i]) == 0:

                            # Determine if walker progress is wihin the threshold.
                            if self.best_prog == 'min':
                                if walker_progress[i] < threshold:
                                    max_tups.append((value, new_walker_weights[i], walker_progress[i], i))
                            elif self.best_prog == 'max':
                                if walker_progress[i] > threshold:
                                    max_tups.append((value, new_walker_weights[i], walker_progress[i], i))

            #print('Walker to be cloned is:', max_tups)

            
            if len(max_tups) > 0:
                max_value, max_weight, progress, max_idx = max(max_tups)                
            pot_merge_pairs = self._find_eligible_merge_pairs(new_walker_weights, distance_matrix, max_idx, new_num_walker_copies)

            merge_pair= self._calc_variation_loss(walker_variations, new_walker_weights, pot_merge_pairs)

            # did we find a closewalk?
            #condition_list = np.array([i is not None for i in [min_idx, max_idx, closewalk]])
            #if we find a walker for cloning, a walker and its close neighbor for merging
            if len(merge_pair) > 0 and len(max_tups) > 0:
                min_idx = merge_pair [0]
                closewalk = merge_pair [1]
                
                # change new_amp
                tempsum = new_walker_weights[min_idx] + new_walker_weights[closewalk]
                new_num_walker_copies[min_idx] = new_walker_weights[min_idx]/tempsum
                new_num_walker_copies[closewalk] = new_walker_weights[closewalk]/tempsum
                new_num_walker_copies[max_idx] += 1

                # re-determine variation function, and walker_variations values
                new_variation, walker_variations = self._calcvariation(new_walker_weights, new_num_walker_copies, distance_matrix)

                if new_variation > variation:
                    variations.append(new_variation)

                    logging.info(f"Variance move to, {new_variation}, accepted")


                    productive = True
                    variation = new_variation

                    # make a decision on which walker to keep
                    # (min_idx, or closewalk), equivalent to:
                    # `random.choices([closewalk, min_idx],
                    #                 weights=[new_walker_weights[closewalk], new_walker_weights[min_idx])`
                    r = rand.uniform(0.0, new_walker_weights[closewalk] + new_walker_weights[min_idx])

                     # keeps closewalk and gets rid of min_idx
                    if r < new_walker_weights[closewalk]:
                        keep_idx = closewalk
                        squash_idx = min_idx

                    # keep min_idx, get rid of closewalk
                    else:
                        keep_idx = min_idx
                        squash_idx = closewalk

                    # update weight
                    new_walker_weights[keep_idx] += new_walker_weights[squash_idx]
                    new_walker_weights[squash_idx] = 0.0

                    # update new_num_walker_copies
                    new_num_walker_copies[squash_idx] = 0
                    new_num_walker_copies[keep_idx] = 1

                    # add the squash index to the merge group
                    merge_groups[keep_idx].append(squash_idx)

                    # add the indices of the walkers that were already
                    # in the merge group that was just squashed
                    merge_groups[keep_idx].extend(merge_groups[squash_idx])

                    # reset the merge group that was just squashed to empty
                    merge_groups[squash_idx] = []

                    # increase the number of clones that the cloned
                    # walker has
                    walker_clone_nums[max_idx] += 1

                    # new variation for starting new stage
                    new_variation, walker_variations = self._calcvariation(new_walker_weights,
                                                                          new_num_walker_copies,
                                                                          distance_matrix)
                    variations.append(new_variation)

                    logging.info(f"variance after selection: {new_variation}")


                # if not productive
                else:
                    new_num_walker_copies[min_idx] = 1
                    new_num_walker_copies[closewalk] = 1
                    new_num_walker_copies[max_idx] -= 1

        # given we know what we want to clone to specific slots
        # (squashing other walkers) we need to determine where these
        # squashed walkers will be merged
        walker_actions = self.assign_clones(merge_groups, walker_clone_nums)

        # because there is only one step in resampling here we just
        # add another field for the step as 0 and add the walker index
        # to its record as well
        for walker_idx, walker_record in enumerate(walker_actions):
            walker_record['step_idx'] = np.array([0])
            walker_record['walker_idx'] = np.array([walker_idx])

        return walker_actions, variations[-1]


    def resample(self, walkers):
        """Resamples walkers based on REVO algorithm
        
        Parameters
        ----------
        
        walkers : list of walkers
        
        Returns
        -------
        
        resampled_walkers : list of resampled_walkers
        
        resampling_data : list of dict of str: value
            The resampling records resulting from the decisions.
        
        resampler_data :list of dict of str: value
            The resampler records resulting from the resampler actions.
        """

        #initialize the parameters
        num_walkers = len(walkers)
        walker_weights = [walker.weight for walker in walkers]
        num_walker_copies = [1 for i in range(num_walkers)]

        # calculate distance matrix
        distance_matrix, images = self._all_to_all_distance(walkers)

        # Calculate the walker progress
        walker_prog = self._calc_progress(walkers)

        logging.info("distance_matrix")
        logging.info("\n{}".format(str(np.array(distance_matrix))))

        # determine cloning and merging actions to be performed, by
        # maximizing the variation, i.e. the Decider
        resampling_data, variation = self.decide(walker_weights, num_walker_copies, distance_matrix, walker_prog)

        # convert the target idxs and decision_id to feature vector arrays
        for record in resampling_data:
            record['target_idxs'] = np.array(record['target_idxs'])
            record['decision_id'] = np.array([record['decision_id']])

        # actually do the cloning and merging of the walkers
        resampled_walkers = self.DECISION.action(walkers, [resampling_data])

       # flatten the distance matrix and give the number of walkers
        # as well for the resampler data, there is just one per cycle
        resampler_data = [{'distance_matrix' : np.ravel(np.array(distance_matrix)),
                           'num_walkers' : np.array([len(walkers)]),
                           'variation' : np.array([variation]),
                         }]

        return resampled_walkers, resampling_data, resampler_data
