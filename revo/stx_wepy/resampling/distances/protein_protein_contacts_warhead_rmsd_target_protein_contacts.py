import numpy as np

from wepy.resampling.distances.distance import Distance
from stx_wepy.resampling.distances.protein_protein_contacts import ProteinProteinContacts
from stx_wepy.resampling.distances.rebinding import RebindingDistance

class ProteinProteinContactsWarheadRMSDProteinProtacContacts(Distance):
    """Distance metric for measuring differences between walker states in                                                        
    regards to the number of protein-protein contacts and the warhead RMSD.                                                                                        
                                                                                                                                 
    Images are produced using the ReceptorDistance.image method. The                            
    distance between images then is the relative difference between                                
    the number of contacts.

    Strength of the contact between the residues between the two protein
    are determined using a sigmoidal function of form y = -1 / (1 + exp(k*(x - x0)))

    where:
    y is the strength of the residue contact.
    x is the minimum distance between the two residues.
    x0 is the distance cutoff (where we want the strength to be at 50%.
    k is the steepness of the curve.
                                                                                     
    """

    def __init__(self,
                 ref_state = None,
                 ligand_idxs = None,
                 binding_site_idxs = None,
                 native_ligand_idxs = None,
                 native_binding_site_idxs = None,
                 prot_1_resids = None,
                 prot_2_resids = None,
                 target_resids = None,
                 protac_resids = None,
                 prot_1_idxs = None,
                 prot_2_idxs = None,
                 target_idxs = None,
                 protac_idxs = None,
                 trajectory = None,
                 distance_cutoff_1 = 0.5,
                 distance_cutoff_2 = 0.5,
                 k_1 = 1,
                 k_2 = 1,
                 warhead_rmsd_weight = 1,
                 protein_protein_contact_dist_weight = 1,
                 target_protact_contacts_dist_weight = 1,
                 **kwargs):


        self.warhead_distance = RebindingDistance(ligand_idxs = ligand_idxs,
                                                  binding_site_idxs = binding_site_idxs,
                                                  native_ligand_idxs = native_ligand_idxs,
                                                  native_binding_site_idxs = native_binding_site_idxs,
                                                  ref_state = ref_state)

        self.protein_protein_contact_distance = ProteinProteinContacts(prot_1_resids = prot_1_resids,
                                                                       prot_2_resids = prot_2_resids,
                                                                       prot_1_idxs = prot_1_idxs,
                                                                       prot_2_idxs = prot_2_idxs,
                                                                       trajectory = trajectory,
                                                                       distance_cutoff = distance_cutoff_1,
                                                                       k = k_1,
                                                                       native_state = ref_state)

        print(target_idxs)
        print(protac_idxs)
        
        self.target_protac_contact_distance = ProteinProteinContacts(prot_1_resids = target_resids,                                                                                                  
                                                                     prot_2_resids = protac_resids,                                                                                                     
                                                                     prot_1_idxs = target_idxs,                                                                                                         
                                                                     prot_2_idxs = protac_idxs,                                                                                                         
                                                                     trajectory = trajectory,                                                                                                           
                                                                     distance_cutoff = distance_cutoff_2,                                                                                               
                                                                     k = k_2,                                                                                                                           
                                                                     native_state = ref_state)


        self.warhead_rmsd_weight = warhead_rmsd_weight
        self.protein_protein_contact_dist_weight = protein_protein_contact_dist_weight
        self.target_protact_contacts_dist_weight = target_protact_contacts_dist_weight



    def image(self, state):
        """ Transform the state into the two images formed from the rebinding distance metric
        and the ProteinProein contacts metric.
        
        Parameters
        ----------
        state : object implementing WalkerState
            State with 'positions' (Nx3 dims) and 'box_vectors' (3x3 array)
            attributes.
        Returns
        -------
        receptor_image : array of float
            The positions of binding site and ligand after
            preprocessing.
        """

        warhead_image = self.warhead_distance.image(state)

        proten_protein_contact_image = self.protein_protein_contact_distance.image(state)

        target_protac_contact_image = self.target_protac_contact_distance.image(state)

        images = np.array([warhead_image, proten_protein_contact_image, target_protac_contact_image])


        
        return images


        
    
    def image_distance(self, image_a, image_b):

        # Seperate images
        warhead_image_a = image_a[0]
        protein_protein_contact_image_a = image_a[1]
        target_protac_contact_image_a = image_a[2]

        warhead_image_b = image_b[0]
        protein_protein_contact_image_b = image_b[1]
        target_protac_contact_image_b = image_b[2]

        print(target_protac_contact_image_a)
        print(target_protac_contact_image_b)



        # warhead image distance
        warhead_distance = self.warhead_distance.image_distance(warhead_image_a,
                                                                warhead_image_b)


        # warhead image distance                                                                                 
        protein_protein_contact_distance = self.protein_protein_contact_distance.image_distance(protein_protein_contact_image_a,
                                                                                                protein_protein_contact_image_b)

        protein_protein_contact_distance = self.target_protac_contact_distance.image_distance(target_protac_contact_image_a,                                                                          
                                                                                              target_protac_contact_image_b)                                                                           



        distance = (self.warhead_rmsd_weight * warhead_distance +
                    self.protein_protein_contact_dist_weight * protein_protein_contact_distance +
                    self.target_protact_contacts_dist_weight * protein_protein_contact_distance)
        
        return distance
