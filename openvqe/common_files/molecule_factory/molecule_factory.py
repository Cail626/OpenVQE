import numpy as np
import scipy
from numpy import binary_repr
from qat.fermion import ElectronicStructureHamiltonian
from qat.fermion.chemistry.pyscf_tools import perform_pyscf_computation
from qat.fermion.chemistry.ucc import (
    convert_to_h_integrals,
    transform_integrals_to_new_basis,
)
from qat.fermion.chemistry.ucc_deprecated import (
    get_active_space_hamiltonian,
    get_cluster_ops_and_init_guess,
)
from qat.fermion.transforms import (
    get_bk_code,
    get_jw_code,
    get_parity_code,
    recode_integer,
    transform_to_bk_basis,
    transform_to_jw_basis,
    transform_to_parity_basis,
)

from .generator.generator_excitations import (
    singlet_gsd,
    singlet_sd,
    singlet_upccgsd,
    spin_complement_gsd,
    spin_complement_gsd_twin,
    uccsd,
)


class MoleculeFactory:

    def __init__(
        self, molecule, active, transform='JW', threshold_1=None, threshold_2=None
    ):
        """Generate a molecule

        Args:
            molecule (Molecule): molecule to study
            active (bool): active space selection
            transform (str, optional): _description_. Defaults to 'JW'.
            threshold_1 (_type_, optional): _description_. Defaults to None.
            threshold_2 (_type_, optional): _description_. Defaults to None.
        """
        self.molecule = molecule
        self.active = active
        self.transform=transform
        self.threshold_1=threshold_1
        self.threshold_2=threshold_2
        
    def build(self, display=True):
        hamiltonian, hamiltonian_sp, n_elec, noons_full, orb_energies_full, info = self.generate_hamiltonian(display=True)
        hf_init = self.find_hf_init(hamiltonian, n_elec, noons_full, orb_energies_full)
        reference_ket, hf_init_sp = self.get_reference_ket(hf_init, hamiltonian_sp.nbqbits, self.transform)
        self.hf_init_sp = hf_init_sp
        self.hamiltonian_sp = hamiltonian_sp
        
    def generate_hamiltonian(self, display=True):
        (
            rdm1,
            orbital_energies,
            nuclear_repulsion,
            n_elec,
            one_body_integrals,
            two_body_integrals,
            info,
        ) = perform_pyscf_computation(
            geometry=self.geometry, basis=self.basis, spin=self.spin, charge=self.charge, run_fci=True
        )
        if display:
            print("Number of electrons = ", n_elec)
            print(
                "Number of qubits before active space selection = ", rdm1.shape[0] * 2
            )
            # print("rdm1", rdm1)
            # print(info)
            print("Orbital energies = ", orbital_energies)
            print("Nuclear repulsion = ", nuclear_repulsion)

        # nuclear_repulsion = molbasis.nuclear_repulsion

        hpq, hpqrs = convert_to_h_integrals(one_body_integrals, two_body_integrals)
        if self.active == False:

            hamiltonian = ElectronicStructureHamiltonian(
                hpq, hpqrs, constant_coeff=nuclear_repulsion
            )
            noons, basis_change = np.linalg.eigh(rdm1)
            noons = list(reversed(noons))
            if display:
                print("Noons = ", noons)
            noons_full, orb_energies_full = [], []
            for ind in range(len(noons)):
                noons_full.extend([noons[ind], noons[ind]])
                orb_energies_full.extend([orbital_energies[ind], orbital_energies[ind]])
            hamiltonian_sp = None
            if transform == "JW":
                trafo, code = transform_to_jw_basis, get_jw_code
                hamiltonian_sp = trafo(hamiltonian)
            elif transform == "Bravyi-Kitaev":
                trafo, code = transform_to_bk_basis, get_bk_code
                hamiltonian_sp = trafo(hamiltonian)
            elif transform == "parity_basis":
                trafo, code = transform_to_parity_basis, get_parity_code
                hamiltonian_sp = trafo(hamiltonian)
            return (
                hamiltonian,
                hamiltonian_sp,
                n_elec,
                noons_full,
                orb_energies_full,
                info,
            )

        noons, basis_change = np.linalg.eigh(rdm1)
        noons = list(reversed(noons))
        if display:
            print("Noons = ", noons)
        # noons, basis_change = np.linalg.eigh(rdm1)
        # noons = list(reversed(noons))  # need to put noons in decreasing order
        basis_change = np.flip(basis_change, axis=1)
        one_body_integrals, two_body_integrals = transform_integrals_to_new_basis(
            one_body_integrals, two_body_integrals, basis_change
        )

        #         threshold_1 = 2-(noons[0]+noons[1])/2
        threshold_1 = 2 - noons[0]
        if len(noons) < 3:
            threshold_2 = 0.01
        else:
            threshold_2 = noons[3]
        if display:
            print("threshold_1 chosen = ", threshold_1)
            print("threshold_2 chosen = ", threshold_2)
        hamiltonian_active, active_inds, occ_inds = get_active_space_hamiltonian(
            one_body_integrals,
            two_body_integrals,
            noons,
            n_elec,
            nuclear_repulsion,
            threshold_1=threshold_1,
            threshold_2=threshold_2,
        )
        if display:
            print(
                "Number of qubits after active space selection =",
                hamiltonian_active.nbqbits,
            )
        # print(hamiltonian_active)
        # print(geek.identity(4))
        # hamiltonian_active = hamiltonian_active -(6.896385617948947*geek.identity(4))
        ## print('hamiltonian_active:', hamiltonian_active)
        active_noons, active_orb_energies = [], []
        for ind in active_inds:
            active_noons.extend([noons[ind], noons[ind]])
            active_orb_energies.extend([orbital_energies[ind], orbital_energies[ind]])
        nb_active_els = n_elec - 2 * len(occ_inds)
        if display:
            print("length of active noons: ", len(active_noons))
            print("length of orbital energies: ", len(active_orb_energies))
        # import numpy as np
        # eigen,_ = np.linalg.eigh(hamiltonian_active.get_matrix())
        # print("eigen",eigen)

        # print("eigen",min(eigen))
        hamiltonian_active_sp = None
        if transform == "JW":
            trafo, code = transform_to_jw_basis, get_jw_code
            hamiltonian_active_sp = trafo(hamiltonian_active)
        elif transform == "Bravyi-Kitaev":
            trafo, code = transform_to_bk_basis, get_bk_code
            hamiltonian_active_sp = trafo(hamiltonian_active)
        elif transform == "parity_basis":
            trafo, code = transform_to_parity_basis, get_parity_code
            hamiltonian_active_sp = trafo(hamiltonian_active)
        return (
            hamiltonian_active,
            hamiltonian_active_sp,
            nb_active_els,
            active_noons,
            active_orb_energies,
            info,
        )

    def calculate_uccsd(self, molecule_symbol, transform, active):
        if active == False:
            (
                hamiltonian,
                hamiltonian_sp,
                n_elec,
                noons_full,
                orb_energies_full,
                info,
            ) = self.generate_hamiltonian(
                molecule_symbol, active=False, transform=transform, display=False
            )
            pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init = uccsd(
                hamiltonian, n_elec, noons_full, orb_energies_full, transform
            )
            return pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init
        else:
            (
                hamiltonian_active,
                hamiltonian_active_sp,
                nb_active_els,
                active_noons,
                active_orb_energies,
                info,
            ) = self.generate_hamiltonian(
                molecule_symbol, active=True, transform=transform
            )
            pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init = uccsd(
                hamiltonian_active,
                nb_active_els,
                active_noons,
                active_orb_energies,
                transform,
            )
            return pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init

    def find_hf_init(self, hamiltonian, n_elec, noons_full, orb_energies_full):
        _, _, hf_init = get_cluster_ops_and_init_guess(
            n_elec, noons_full, orb_energies_full, hamiltonian.hpqrs
        )
        return hf_init

    # TESTING HOW TO GET THE REFERENCE KET
    def get_reference_ket(self, hf_init, nbqbits):
        if self.transform == "JW":
            hf_init_sp = recode_integer(hf_init, get_jw_code(nbqbits))
        elif self.transform == "Bravyi-Kitaev":
            hf_init_sp = recode_integer(hf_init, get_bk_code(nbqbits))
        elif self.transform == "parity_basis":
            hf_init_sp = recode_integer(hf_init, get_parity_code(nbqbits))
        ket_hf = binary_repr(hf_init_sp)
        list_ket_hf = [int(c) for c in ket_hf]
        reference_state = self.from_ket_to_vector(list_ket_hf)
        sparse_reference_state = scipy.sparse.csr_matrix(
            reference_state, dtype=complex
        ).transpose()
        return sparse_reference_state, hf_init_sp

    def from_ket_to_vector(self, ket):
        state_vector = [1]
        for i in ket:
            qubit_vector = [not i, i]
            state_vector = np.kron(state_vector, qubit_vector)
        return state_vector

    def generate_cluster_ops(
        self, molecule_symbol, type_of_generator, transform, active=False
    ):
        r, geometry, charge, spin, basis = self.get_parameters(molecule_symbol)
        (
            rdm1,
            orbital_energies,
            nuclear_repulsion,
            n_elec,
            one_body_integrals,
            two_body_integrals,
            info,
        ) = perform_pyscf_computation(
            geometry=geometry, basis=basis, spin=spin, charge=charge, run_fci=True
        )

        orbital_number = len(orbital_energies)

        if active == True:
            (
                hamiltonian_active,
                hamiltonian_active_sp,
                nb_active_els,
                active_noons,
                active_orb_energies,
                info,
            ) = self.generate_hamiltonian(
                molecule_symbol, active=True, transform=transform, display=False
            )
            # print('here: ', len(active_orb_energies))
            orbital_number = int(len(active_orb_energies) / 2)
            n_elec = nb_active_els

        # orb_energies_full = orbital_energies
        pool_size, cluster_ops, cluster_ops_sp = None, None, None
        if type_of_generator == "singlet_sd":
            pool_size, cluster_ops, cluster_ops_sp = singlet_sd(
                n_elec, orbital_number, transform
            )
        elif type_of_generator == "singlet_gsd":
            pool_size, cluster_ops, cluster_ops_sp = singlet_gsd(
                n_elec, orbital_number, transform
            )
        elif type_of_generator == "spin_complement_gsd":
            pool_size, cluster_ops, cluster_ops_sp = spin_complement_gsd(
                n_elec, orbital_number, transform
            )
        elif type_of_generator == "spin_complement_gsd_twin":
            pool_size, cluster_ops, cluster_ops_sp = spin_complement_gsd_twin(
                n_elec, orbital_number, transform
            )
        elif type_of_generator == "sUPCCGSD":
            # user can type here perm = 1,2,3,4,....10
            perm = 2
            pool_size, cluster_ops, cluster_ops_sp = singlet_upccgsd(
                orbital_number, transform,perm
            )
        elif type_of_generator == "QUCCSD":
            (
                pool_size,
                cluster_ops,
                cluster_ops_sp,
                theta_MP2,
                hf_init,
            ) = self.calculate_uccsd(molecule_symbol, transform, active=active)
            return pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init
        elif type_of_generator == "UCCSD":
            (
                pool_size,
                cluster_ops,
                cluster_ops_sp,
                theta_MP2,
                hf_init,
            ) = self.calculate_uccsd(molecule_symbol, transform, active=active)
            return pool_size, cluster_ops, cluster_ops_sp, theta_MP2, hf_init
        else:
            return
        return pool_size, cluster_ops, cluster_ops_sp
