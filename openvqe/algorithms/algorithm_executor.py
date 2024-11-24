from .adapt.qubit_adapt_vqe import qubit_adapt_vqe
from .ucc_family.get_energy_ucc import EnergyUCC
from .ucc_family.get_energy_qucc import EnergyQUCC

class AlgorithmExecutor:
     
    def __init__(name):
        supported_algorithm = {
            "UCC": EnergyUCC().get_energies,
            "QUCC": EnergyQUCC().get_energies,
            "ADAPT-VQE": qubit_adapt_vqe
        }
        
        algorithm = supported_algorithm[name]
        if not algorithm:
            raise ValueError(f'algorithm {name} is not supported. The currently supported algorithms are: {supported_algorithm.keys()}')
        
        return algorithm