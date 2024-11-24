import numpy as np

def get_parameters(self, molecule_symbol):
    """This method will be used to multiply two numbers
    :param string molecule_symbol: The symbol of the molecule
    :returns: r, geometry, charge, spin, basis
    :rtype: multiple
    """
    if molecule_symbol == "LIH":
        r = 1.45
        geometry = [("Li", (0, 0, 0)), ("H", (0, 0, r))]
        charge = 0
        spin = 0
        basis = "sto-3g"   
    elif molecule_symbol == "H2":
        r = 0.75
        geometry = [("H", (0, 0, 0)), ("H", (0, 0, r))]
        charge = 0
        spin = 0
        basis = "6-31g"
    elif molecule_symbol == "H4":
        # H4
        r = 0.85
        geometry = [
            ("H", (0, 0, 0)),
            ("H", (0, 0, 1 * r)),
            ("H", (0, 0, 2 * r)),
            ("H", (0, 0, 3 * r)),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "H6":
        r = 1.0
        geometry = [
            ("H", (0, 0, 0)),
            ("H", (0, 0, 1 * r)),
            ("H", (0, 0, 2 * r)),
            ("H", (0, 0, 3 * r)),
            ("H", (0, 0, 4 * r)),
            ("H", (0, 0, 5 * r)),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "H8":
        r = 1.0
        geometry = [
            ("H", (0, 0, 0)),
            ("H", (0, 0, 1 * r)),
            ("H", (0, 0, 2 * r)),
            ("H", (0, 0, 3 * r)),
            ("H", (0, 0, 4 * r)),
            ("H", (0, 0, 5 * r)),
            ("H", (0, 0, 6 * r)),
            ("H", (0, 0, 7 * r)),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "H10":
        r = 1.0
        geometry = [
            ("H", (0, 0, 0)),
            ("H", (0, 0, 1 * r)),
            ("H", (0, 0, 2 * r)),
            ("H", (0, 0, 3 * r)),
            ("H", (0, 0, 4 * r)),
            ("H", (0, 0, 5 * r)),
            ("H", (0, 0, 6 * r)),
            ("H", (0, 0, 7 * r)),
            ("H", (0, 0, 8 * r)),
            ("H", (0, 0, 9 * r)),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "BeH2":
        r = 1.4
        geometry = [("Be", (0, 0, 0 * r)), ("H", (0, 0, r)), ("H", (0, 0, -r))]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "HeH+":
        r = 1.0
        geometry = [("He", (0, 0, 0)), ("H", (0, 0, r))]
        charge = 1
        spin = 0
        basis = "6-31g"
    elif molecule_symbol == "HF":
        r = 1.0
        geometry = [("F", (0, 0, 0 * r)), ("H", (0, 0, r))]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "HO":
        r = 1.8
        geometry = [("H", (0, 0, 0 * r)), ("O", (0, 0, 1 * r))]
        charge = -1
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "H2O":
        r = 1.0285
        theta = 0.538 * np.pi
        geometry = [
            ("O", (0, 0, 0 * r)),
            ("H", (0, 0, r)),
            ("H", (0, r * np.sin(np.pi - theta), r * np.cos(np.pi - theta))),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "NH3":
        r = 1.0703
        theta = (100.107 / 180) * np.pi
        geometry = [
            ("N", (0, 0, 0 * r)),
            (
                "H",
                (
                    0,
                    2 * (np.sin(theta / 2) / np.sqrt(3)) * r,
                    np.sqrt(1 - 4 * np.sin(theta / 2) ** 2 / 3) * r,
                ),
            ),
            (
                "H",
                (
                    np.sin(theta / 2) * r,
                    -np.sin(theta / 2) / np.sqrt(3) * r,
                    np.sqrt(1 - 4 * np.sin(theta / 2) ** 2 / 3) * r,
                ),
            ),
            (
                "H",
                (
                    -np.sin(theta / 2) * r,
                    -np.sin(theta / 2) / np.sqrt(3) * r,
                    np.sqrt(1 - 4 * np.sin(theta / 2) ** 2 / 3) * r,
                ),
            ),
        ]
        charge = 0
        spin = 0
        basis = "sto-3g"
    elif molecule_symbol == "CO2":
        r = 1.22
        geometry = [
            ["C", [0.0, 0.0, 8.261342997000753e-07]],
            [
                "O",
                [
                    1.0990287608769004e-18,
                    2.7114450405987004e-19,
                    1.2236575813458745,
                ],
            ],
            [
                "O",
                [
                    2.696319376811295e-22,
                    2.4247676462727696e-23,
                    -1.2236561920609494,
                ],
            ],
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "SO2":
        r = 1.0
        geometry = [
            ("S", (0.0, 0.0, 0.0)),
            ("O", (0.0, 1.2371, 0.7215)),
            ("O", (0.0, -1.2371, 0.7215)),
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "Cl2":
        r = 1.0
        geometry = [("cl", (0.0, 0.0, 0.0)), ("cl", (0.0, 0.0, 1.9879))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "S2":
        r = 1.0
        geometry = [("S", (0.0, 0.0, 0.0)), ("S", (0.0, 0.0, 1.8892))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "C2H2":
        r = 1.0
        geometry = [
            ("C", (0.0, 0.0, 0.6063)),
            ("C", (0.0, 0.0, -0.6063)),
            ("H", (0.0, 0.0, 1.6941)),
            ("H", (0.0, 0.0, -1.6941)),
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "CO":
        r = 1.0
        geometry = [("C", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.1282))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "N2":
        r = 1.0
        geometry = [("N", (0.0, 0.0, 0.5488)), ("N", (0.0, 0.0, -0.5488))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "F2":
        r = 1.0
        geometry = [("F", (0.0, 0.0, 0.0)), ("F", (0.0, 0.0, 1.4119))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "CH4":
        r = 1.0
        geometry = [
            ("C", (0.0, 0.0, 0.0)),
            ("H", (0.6276, 0.6276, 0.6276)),
            ("H", (0.6276, -0.6276, -0.6276)),
            ("H", (-0.6276, 0.6276, -0.6276)),
            ("H", (-0.6276, -0.6276, 0.6276)),
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "C2H4":
        r = 1.0
        geometry = [
            ("C", (0.0, 0.0, 0.6695)),
            ("C", (0.0, 0.0, -0.6695)),
            ("H", (0.0, 0.9289, 1.2321)),
            ("H", (0.0, -0.9289, 1.2321)),
            ("H", (0.0, 0.9289, -1.2321)),
            ("H", (0.0, -0.9289, -1.2321)),
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "CHN":
        r = 1.0
        geometry = [
            ("C", (0.0, 0.0, 0.0)),
            ("H", (0.0, 0.0, 1.0640)),
            ("N", (0.0, 0.0, -1.1560)),
        ]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "O2":
        r = 1.0
        geometry = [("O", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.2075))]
        basis = "sto-3g"
        spin = 0
        charge = 0
    elif molecule_symbol == "NO":
        r = 1.0
        geometry = [("N", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.1508))]
        basis = "sto-3g"
        spin = 0
        charge = 1
    return r, geometry, charge, spin, basis