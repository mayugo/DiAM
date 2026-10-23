def f_fatigue_modified_Se(surface, d_eq, load_type, Sut, Sepr):
    """Calculate Se = ka·kb·kc·Se'.

    surface   : surface factor ka (options: 'ground', 'machined', 'hotrolled',
                'forged', 'casted')
    d_eq      : equivalent diameter in mm, for the size factor kb
    load_type : load factor kc (options: 'bending', 'axial', 'torsion')
                NOTE: for 'torsion' kc = 1 is only valid when the fatigue
                analysis is done with the von Mises equivalent stress. For pure
                torsion analysed with shear stresses, Shigley gives kc = 0.59.
    Sut       : ultimate tensile strength in MPa
    Sepr      : Se', value in MPa or material for estimation
                (options: 'steel', 'castiron', 'aluminium', 'copper', 'polymer')

    Returns (Se, ka, kb, kc).
    """
    # --- Surface factor ka = a·Sut^b ---
    surface_coef = {
        'ground':    (1.58, -0.085),
        'machined':  (4.51, -0.265),
        'hotrolled': (57.7, -0.718),
        'forged':    (272.0, -0.995),
        'casted':    (1.0, 0.0),
    }
    if surface not in surface_coef:
        raise ValueError(f"f_fatigue_modified_Se: surface option not valid: "
                         f"'{surface}'. Options: {', '.join(surface_coef)}")
    a, b = surface_coef[surface]
    ka = a * Sut**b

    # --- Size factor kb ---
    if load_type == 'axial':
        kb = 1.0
    elif load_type in ('bending', 'torsion'):
        if d_eq < 2.79:
            kb = (2.79 / 7.62)**(-0.107)
        elif d_eq > 51:
            kb = 1.51 * d_eq**(-0.157)
        else:
            kb = (d_eq / 7.62)**(-0.107)
    else:
        raise ValueError(f"f_fatigue_modified_Se: load_type option not valid: "
                         f"'{load_type}'. Options: bending, axial, torsion")

    # --- Load factor kc ---
    if load_type == 'bending':
        kc = 1.0
    elif load_type == 'axial':
        kc = 0.85
    else:  # 'torsion'
        kc = 1.0   # only valid with von Mises equivalent stress (see NOTE)

    # --- Endurance limit of the test specimen Se' ---
    if isinstance(Sepr, str):
        if Sepr == 'steel':
            Sepr = 0.5 * Sut if Sut < 1400 else 700.0
        elif Sepr == 'castiron':
            Sepr = 0.4 * Sut if Sut < 400 else 160.0
        elif Sepr == 'aluminium':
            Sepr = 0.4 * Sut if Sut < 330 else 130.0
        elif Sepr == 'copper':
            Sepr = 0.4 * Sut if Sut < 280 else 100.0
        elif Sepr == 'polymer':
            Sepr = 0.2 * Sut  # to Sepr = 0.4*Sut
        else:
            raise ValueError(f"f_fatigue_modified_Se: Sepr option not valid: "
                             f"'{Sepr}'. Introduce a value in MPa or one of: "
                             f"steel, castiron, aluminium, copper, polymer")

    Se = Sepr * ka * kb * kc
    return Se, ka, kb, kc
