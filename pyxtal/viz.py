"""
This module handles visualization. Mostly powered by py3Dmol
"""


def addBox(view, vecs, viewer=None):
    pts = [
           ([0,0,0], vecs[0,:]),
           ([0,0,0], vecs[1,:]),
           ([0,0,0], vecs[2,:]),
           (vecs[0,:], vecs[0,:]+vecs[1,:]),
           (vecs[0,:], vecs[0,:]+vecs[2,:]),
           (vecs[1,:], vecs[0,:]+vecs[1,:]),
           (vecs[1,:], vecs[2,:]+vecs[1,:]),
           (vecs[2,:], vecs[0,:]+vecs[2,:]),
           (vecs[2,:], vecs[1,:]+vecs[2,:]),
           (vecs[0,:]+vecs[1,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
           (vecs[0,:]+vecs[2,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
           (vecs[1,:]+vecs[2,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
          ]

    for i, pt in enumerate(pts):
        pt1, pt2 = pt

        if viewer is None:
            view.addLine({"start": {"x":pt1[0], "y": pt1[1], "z": pt1[2]}, 
                           "end":  {"x":pt2[0], "y": pt2[1], "z": pt2[2]}})
        else:
            view.addLine({"start": {"x":pt1[0], "y": pt1[1], "z": pt1[2]}, 
                           "end":  {"x":pt2[0], "y": pt2[1], "z": pt2[2]}},
                           viewer = viewer)
        
    # add unitcell labels
    labels = {"o": [0,0,0],
              "a": vecs[0,:],
              "b": vecs[1,:],
              "c": vecs[2,:]}
    for key in labels.keys():
        text, pos = key, labels[key]
        if viewer is None:
            view.addLabel(text, {"position": {"x": pos[0], "y": pos[1], "z": pos[2]}, 
                                 "fontColor":"black",
                                 "backgroundColor": "white", 
                                 "fontsize": 12,
                                 "backgroundOpacity": "0.1",
                                })
        else:
            view.addLabel(text, {"position": {"x": pos[0], "y": pos[1], "z": pos[2]}, 
                                 "fontColor":"black",
                                 "backgroundColor": "white", 
                                 "fontsize": 12,
                                 "backgroundOpacity": "0.1",
                                },
                              viewer = viewer)
        
def addlines(view, orig, axes, viewer=None):
    for ax in axes:
        ax += orig
        if viewer is None:
            view.addArrow({"start": {"x":orig[0], "y": orig[1], "z":orig[2]}, 
                       "end": {"x":ax[0], "y":ax[1], "z":ax[2]}})
        else:
            view.addArrow({"start": {"x":orig[0], "y": orig[1], "z":orig[2]}, 
                           "end": {"x":ax[0], "y":ax[1], "z":ax[2]}},
                           viewer=viewer)



def display_atomic(struc, size=(600,300), scale=0.25, radius=0.10, supercell=(1,1,1), show_wp=False):
    """
    display the molecular crystals generated from pyxtal. If the animation is False, 
    only dump the structure to cif and let py3Dmol display it. If animation is on,
    show the generation of molecules in the crystal as an animation.

    Args:
        size: (width, height) in tuple
        scale: the size of sphere
        radius: the size of stick
        supercell: replicate the crystal (valid only when animation is False)
        show_wp: whether or not highlight the unique wyckoff site

    Returns:
        py3Dmol object
    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)
    if struc.dim == 0:
        fmt = "xyz"
    else:
        fmt = "cif"

    txt = struc.to_file(fmt=fmt)
    view.addModel(txt, fmt, {'doAssembly':True,'duplicateAssemblyAtoms':True})
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':scale},
                   'stick':{'colorscheme':'Jmol', 'radius':radius}})
    if struc.dim != 0:
        view.addUnitCell()
        A, B, C = supercell
        view.replicateUnitCell(A,B,C)
        if show_wp:
            view.setStyle({'sym':2},{'sphere':{'scale':scale*1.1,'color':'blue'}})
    return view.zoomTo()

def display_molecular(struc, size=(600, 300), supercell=(1,1,1), axis=None, animation=False, interval=2000):
    """
    display the molecular crystals generated from pyxtal. If the animation is False, 
    only dump the structure to cif and let py3Dmol display it. If animation is on,
    show the generation of molecules in the crystal as an animation.

    Args:
        size: (width, height) in tuple
        supercell: replicate the crystal (valid only when animation is False)
        axis: 3-vector to reprent the rotational axis
        animation: whether or not display the animation

    Returns:
        py3Dmol object
    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)
    if not animation:
        cif = struc.to_file()
        view.addModel(cif,'cif',{'doAssembly':True,
                                 'duplicateAssemblyAtoms':True,
                                 'normalizeAssembly':True})
        view.setStyle({'stick':{'colorscheme':'greenCarbon'}})
        view.addUnitCell()
        if axis is not None:
            site = struc.mol_sites[0]
            for id in range(len(site.wp)):
                mol = site.get_mol_object(id)
                addlines(view, mol.cart_coords.mean(axis=0), [axis.copy()])
 
        A, B, C = supercell
        view.replicateUnitCell(A,B,C)
    else:
        cifs = ""
        for i in range(len(struc.mol_sites[0].wp)):
            cifs += struc.to_file(sym_num=i+1)
        view.addModelsAsFrames(cifs, 'cif')
        view.animate({'loop': 'forward', 'interval': interval})
        view.addUnitCell()
        view.setStyle({'model':0},{'stick':{'colorscheme':'greenCarbon'}})
    return view.zoomTo()

def display_molecular_site(site, id=None, size=(400, 300), axis=True, ax_id=range(3), box=False):
    """
    display the Wyckoff site in the molecular crystals generated from pyxtal. 

    Args:
        site: pyxtal.wyckoff_site.mol_site object
        id: list of molecules to display. If None, display all molecules in this site
        size: (width, height) in tuple
        axis: boolean, whether or not display the rotational axis
        axis_label: whether or not display the unitcell axis labels

    Returns:
        py3Dmol object
    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)

    if id is None:
        ids = range(site.wp.multiplicity)
    else:
        ids = [id]

    for i in ids:
        mol = site.get_mol_object(i)
        cell, _, center = site.molecule.get_box_coordinates(mol.cart_coords)
        view.addModel(mol.to(fmt='xyz'), 'xyz')
        if axis:
            axes = site.molecule.get_principle_axes(mol.cart_coords)
            addlines(view, mol.cart_coords.mean(axis=0), axes.T[ax_id]*5, viewer=(0,1))
        if box:
            center_spec = {'x':center[0], 'y':center[1], 'z': center[2]}
            w_spec = {'x':cell[0,0], 'y': cell[0,1], 'z': cell[0,2]}
            h_spec = {'x':cell[1,0], 'y': cell[1,1], 'z': cell[1,2]}
            d_spec = {'x':cell[2,0], 'y': cell[2,1], 'z': cell[2,2]}
            view.addBox({'center': center_spec,
                        'dimensions': {'w': w_spec, 'h': h_spec, 'd': d_spec},
                        'color':'magenta',
                        'alpha': 0.5,
                       })
   
    addBox(view, site.lattice.matrix, viewer=(0,1))
    view.setStyle({'stick':{'colorscheme':'greenCarbon'}})

    return view.zoomTo()

def display_molecules(molecules, size=(400,300), animation=False, center=False):
    """
    display the molecules in Pymatgen object.

    Args:
        molecules: a list of pymatgen molecules
        size: (width, height) in tuple
        animation: whether or not display animation
        center: highlight center or not

    Returns:
        py3Dmol object

    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)
    mol_strs = ""
    for mol in molecules:
        mol_strs += mol.to(fmt='xyz') + '\n'
    if animation:
        view.addModelsasFrames(mol_strs, 'xyz')
        view.animate({'loop': 'forward'})
    else:
        view.addModels(mol_strs, 'xyz')
    view.setStyle({'stick':{'colorscheme':'greenCarbon'}})

    return view.zoomTo()

def display_molecule(molecule, center, cell, size=(400,300)):
    """
    display the molecules in Pymatgen object.

    Args:
        mol: pymatgen molecule
        center: molecular
        size: (width, height) in tuple

    Returns:
        py3Dmol object

    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)
    mol_strs = ""
    mol_strs += molecule.to(fmt='xyz') + '\n'
    view.addModels(mol_strs, 'xyz')
    center_spec = {'x':center[0], 'y':center[1], 'z': center[2]}
    w_spec = {'x':cell[0,0], 'y': cell[0,1], 'z': cell[0,2]}
    h_spec = {'x':cell[1,0], 'y': cell[1,1], 'z': cell[1,2]}
    d_spec = {'x':cell[2,0], 'y': cell[2,1], 'z': cell[2,2]}
    view.addBox({'center': center_spec,
                 'dimensions': {'w': w_spec, 'h': h_spec, 'd': d_spec},
                 'color':'magenta',
                 'opacity': 0.4,
                 #'alpha': 0.6,
                })
    view.setStyle({'stick':{'colorscheme':'greenCarbon'}})

    return view.zoomTo()


def display_mol_crystals(strucs, size=(600, 300), supercell=(1,1,1), axis=None, animation='slider', interval=2000):
    """
    display the molecular crystals generated from pyxtal. 
    two modes of animations are supported: slider or movie

    Args:
        size: (width, height) in tuple
        supercell: replicate the crystal (valid only when animation is False)
        animation: slider or movie

    Returns:
        py3Dmol object
    """
    import py3Dmol

    (width, height) = size
    view = py3Dmol.view(height=height, width=width)
    if animation=='slider':
        from ipywidgets import interact, IntSlider
        def conf_viewer(idx):
            return strucs[idx].show(size=size, supercell=supercell, axis=axis)
        interact(conf_viewer, idx=IntSlider(min=0, max=len(strucs)-1, description='id:'))

    elif animation=='movie':
        cifs = ""
        for struc in strucs:
            cifs += struc.to_file()
        view.addModelsAsFrames(cifs, 'cif')
        view.animate({'loop': 'forward', 'interval': interval})
        view.addUnitCell()
        view.setStyle({'model':0},{'stick':{'colorscheme':'greenCarbon'}})
        return view.zoomTo()
    else:
        raise ValueError("only movie or slider is supported")
       


def display_cluster(molecules, cell, Ps, N_cut=12, size=(400,300), style='sphere'):
    import py3Dmol
    models = {}
    
    view = py3Dmol.view()
    mol_str = ""
    mol_str += molecules[0].to(fmt='xyz') + '\n'
    view.addModel(mol_str, 'xyz')
    model = view.getModel()
    model.setStyle({}, {"sphere": {'colorscheme':'grayCarbon', 
                                   'scale':0.7}})

    center = molecules[0].center_of_mass
    center_spec = {'x':center[0], 'y':center[1], 'z': center[2]}
    w_spec = {'x':cell[0,0], 'y': cell[0,1], 'z': cell[0,2]}
    h_spec = {'x':cell[1,0], 'y': cell[1,1], 'z': cell[1,2]}
    d_spec = {'x':cell[2,0], 'y': cell[2,1], 'z': cell[2,2]}
    #view.addBox({'center': center_spec,
    #             'dimensions': {'w': w_spec, 'h': h_spec, 'd': d_spec},
    #             'color':'magenta',
    #             'alpha': 0.5,
    #             #'opacity': 0.1,
    #            })
    addBox(view, cell)
       
    for i in range(1, len(molecules)):
        mol = molecules[i]
        mol_strs = ""
        mol_strs += mol.to(fmt='xyz') + '\n'
        view.addModel(mol_strs, 'xyz')
        model = view.getModel()
        if sum(Ps) > 0:
            if Ps[i-1] == 0:
                ctype = 1
            else:
                ctype = 2
        else:
            if i > N_cut:
                ctype = 1
            else:
                ctype = 2
        if ctype == 1:
            color = 'greenCarbon'
            opacity = 0.65
        else:
            color = 'cyanCarbon'
            opacity = 0.65
        if style == 'sphere':
            model.setStyle({}, {"sphere": {'colorscheme':color,
                                           'scale':0.5,
                                           'opacity': opacity}})
        else:
            model.setStyle({}, {"stick": {'colorscheme':color, 
                                          'radius': 0.1,
                                          'opacity': opacity}})
    return view.zoomTo({"model": list(models.values())})
