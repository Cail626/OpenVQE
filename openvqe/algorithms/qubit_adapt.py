from openvqe.common_files.qubit_pool import QubitPool
from openvqe.common_files.molecule_factory_with_sparse import MoleculeFactory
from openvqe.adapt.qubit_adapt_vqe import qubit_adapt_vqe

def presentation(molecule_symbol, type_of_generator, transform, active):
    molecule_factory = MoleculeFactory()

    r, geometry, charge, spin, basis = molecule_factory.get_parameters(molecule_symbol)
    print(" --------------------------------------------------------------------------")
    if active:
        print("Running in the active case: ")
    else:
        print("Running in the non active case: ")
    print("                     molecule symbol: %s " %(molecule_symbol))
    print("                     molecule basis: %s " %(basis))
    print("                     type of generator: %s " %(type_of_generator))
    print("                     transform: %s " %(transform))
    print(" --------------------------------------------------------------------------")


def generate_sparse_hamiltonian(molecule_symbol, type_of_generator, transform, active):
    molecule_factory = MoleculeFactory()
    
    print(" --------------------------------------------------------------------------")
    print("                                                          ")
    print("                      Generate Hamiltonians and Properties from :")
    print("                                                          ")
    print(" --------------------------------------------------------------------------")
    print("                                                          ")

    return molecule_factory.generate_hamiltonian(molecule_symbol, active=active, transform=transform)   

def generate_cluster_ops(molecule_symbol, type_of_generator, transform, active):
    molecule_factory = MoleculeFactory()
    print(" --------------------------------------------------------------------------")
    print("                                                          ")
    print("                      Generate Cluster OPS:")
    print("                                                          ")
    print(" --------------------------------------------------------------------------")
    print("                                                           ")

    pool_size,cluster_ops, cluster_ops_sp, cluster_ops_sparse = molecule_factory.generate_cluster_ops(molecule_symbol, type_of_generator=type_of_generator, transform=transform, active=active)

    print('Pool size: ', pool_size)
    print('length of the cluster OP: ', len(cluster_ops))
    print('length of the cluster OPS: ', len(cluster_ops_sp))
    print('length of the cluster _sparse: ', len(cluster_ops_sp))

    return pool_size,cluster_ops,cluster_ops_sp, cluster_ops_sparse

def generate_pool_without_cluster(cluster_ops, nbqbits, molecule_symbol):
    print(" --------------------------------------------------------------------------")
    print("                                                          ")
    print("                      Generate Pool without Cluster:")
    print("                                                          ")
    print(" --------------------------------------------------------------------------")
    print("                                                           ")
    
    qubitpool = QubitPool()
    pool_type = 'random'
    qubit_pool =qubitpool.generate_pool(cluster_ops)
    len_returned_pool, returned_pool = qubitpool.generate_pool_without_cluster(pool_type=pool_type, 
                                                                            nbqbits=nbqbits, 
                                                                            qubit_pool=qubit_pool,
                                                                            molecule_symbol=molecule_symbol)
    return len_returned_pool, returned_pool

def execute(molecule_symbol, type_of_generator, transform, active):
    molecule_factory = MoleculeFactory()

    presentation(molecule_symbol, type_of_generator, transform, active)
    hamiltonian, hamiltonian_sparse, hamiltonian_sp, hamiltonian_sp_sparse, n_elec, noons_full, orb_energies_full, info = generate_sparse_hamiltonian(molecule_symbol, type_of_generator, transform, active)
    pool_size, cluster_ops, cluster_ops_sp, cluster_ops_sparse = generate_cluster_ops(molecule_symbol, type_of_generator, transform, active)
    nbqbits = hamiltonian_sp.nbqbits
    len_returned_pool, returned_pool = generate_pool_without_cluster(cluster_ops, hamiltonian_sp.nbqbits, molecule_symbol)
    hf_init = molecule_factory.find_hf_init(hamiltonian, n_elec, noons_full, orb_energies_full)
    reference_ket, hf_init_sp = molecule_factory.get_reference_ket(hf_init, len(orb_energies_full), transform)
    
    pool_mix = returned_pool
    print("length of the pool",len(pool_mix))
    iterations_sim, iterations_ana, result_sim, result_ana = qubit_adapt_vqe(hamiltonian_sp, hamiltonian_sp_sparse,
    reference_ket, nbqbits, pool_mix, hf_init_sp, info['FCI'],
            adapt_conver    = 'norm',
            adapt_thresh    = 1e-07,
            adapt_maxiter   = 29,
            tolerance_sim = 1e-09,
            method_sim = 'BFGS')
    print("iterations are:",iterations_sim)    
    print("results are:",result_sim)
